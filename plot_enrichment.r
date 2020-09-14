require(tidyr)
require(ggplot2)
require(ggtree)


load("enrichment.Rdata")
load("16S_rooted_tree.Rdata")

#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/plot_enrichment.r --args enrich_file=selected_enrichment.txt

# get absolute path  
#path="/home/jens/CopraRNA-git/coprarna_aux/"
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("plot_enrichment.r","",path)
print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	

# file with modified enrichment
enrich_file="selected_enrichment.txt"


# transforming arguments in valid variables 
args <- commandArgs(trailingOnly = TRUE) 
if(length(args)>0){
	for(i in 1:length(args)){
		temp<-strsplit(args[i],"=")
		temp<-temp[[1]]
		temp1<-temp[1]
		temp2<-temp[2]
		assign(as.character(temp1),temp2)
	}
}

leafs<-function(node, tree){  # finds leafs belonging to an inner node
	topo<-tree[[1]]
	branch<-tree[[2]]
	names(node)<-1
	out<-node
	edge<-c()
	nn<-node
	temp<-c()
	#weight<-c()
	
	while(1>0){
		na<-as.numeric(names(nn))
		n<-which(topo[,1]==nn)
		names(n)<-rep(na,length(n))
		edge<-c(edge,n)
		temp_node<-topo[n,2]
		node<-setdiff(temp_node, node)
		names(node)<-rep(na+1,length(node))
		out<-c(out,node)
		temp<-na.omit(c(temp,node[2:length(node)]))
		nn<-node[1]
		if(is.na(nn)==T & length(temp)>0){
			nn<-temp[1]
			temp<-temp[-1]
		}
		if(is.na(nn)==T & length(temp)==0){
			#print("ff")
			break
			
		}
	}
	dup<-which(duplicated(edge))
	if(length(dup)>0){
		edge<-edge[-dup]
	}
	
	sob<-sum(branch[edge])
	out<-na.omit(out)
	out<-out[which(is.element(out,topo[,1])==F)]
	out
	
}

out2<-read.csv(enrich_file, sep="\t")

term<-as.character(out2[,2])
term2<-as.character(out2[,1])
pval<-as.numeric(out2[,"pVAL"])
node<-setdiff(tree[[1]][,1],tree[[1]][,2])
inner<-unique(tree[[1]][,1])
sel<-which(inner==node)
leaf<-leafs(inner[sel], tree)
leaf<-tree[[3]][leaf]
selection<-leaf

d = fortify(tree)
d = subset(d, isTip)
ord<-with(d, label[order(y, decreasing=T)])
	
nam<-c()
for(i in 1:length(ord)){
	tnam<-grep(gsub("\\..*","",ord[i]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}
nam2<-c()
for(i in 1:length(nam)){
	temp<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp<-paste(temp,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp<-paste(temp, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp)
}
nam2<-paste(nam2, ord, sep="_")


ord2<-c()
for(i in ord){
	p<-grep(i, nam2)
	ord2<-c(ord2,nam2[p])
}




out<-matrix(,0,4)
colnames(out)<-c("term","org","num","pvalue")
for(i in 1:length(term)){
	for(j in selection){
			tmp<-rbind(enrich_list[[j]][[1]][[1]],enrich_list[[j]][[1]][[2]],enrich_list[[j]][[1]][[3]])
			pos<-na.omit(match(term2[i], tmp[,1]))
			tmp2<-c()
			if(length(pos)>0){
				tmp2[1]<-term[i]
				p<-grep(j,nam2)
				tmp2[2]<-nam2[p]
				tmp2[3]<-tmp[pos,"Significant"]
				tmp2[4]<-tmp[pos,"classicFS"]
				out<-rbind(out,tmp2)
			}			
		}
	}
rownames(out)<-1:nrow(out)
	

out<-data.frame(id=as.character(out[,2]),term=as.character(out[,1]), dot_num=as.numeric(out[,3]), pvalue=-1*log10(as.numeric(out[,4])))
	




	p<-as.numeric(out$pvalue)
	breaks<-c(0,1.13,2,3,4,5,6)
	p2<-p
	for(j in breaks){
		tmp<-which(p>=j)
		if(length(tmp)>0){
			p2[tmp]<-j
		
		}
	
	}
	out$colo<-p2
	term_map<-out


col <- c("lightgrey","#7FFFD4","#117893","#34558b","#fd823e","#d13b40" )
col<-c("#117893","#d13b40" )
colfunc <- colorRampPalette(col)
col<-colfunc(6)
col<-c("grey65",col)




p <- ggplot(out, aes((term),(id))) + theme_linedraw()
p + 
 geom_point(aes(size=dot_num, colour=factor(colo))) +
  theme(legend.text.align = 0,axis.text.x = element_text(angle = 90, size=10, vjust=0.2, hjust=1), axis.text.y = element_text(size=6)) +labs(x="", y="", size="term size", colour="p-Value") + 
  scale_y_discrete(limits =rev(ord2)) + scale_x_discrete(limits =(term))  + scale_colour_manual(values = col,labels=c("",expression(""<=5%*%10^{-2}),expression(""<=10^{-2}),expression(""<=10^{-3}),expression(""<=10^{-4}),expression(""<=10^{-5}),expression(""<=10^{-6})))


ggsave("final_enrichment.pdf", height = 3+0.1125*length(tree$tip.label), width = 3.5+0.12*nrow(out2),useDingbats=FALSE)
 
		