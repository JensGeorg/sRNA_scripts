require(tidyr)
require(GOSemSim)
require(org.EcK12.eg.db)
require(ggplot2)
require(ggtree)
require(AnnotationForge)
require(genbankr)
require(topGO)
load("enrichment.Rdata")
load("16S_rooted_tree.Rdata")
weight<-read.csv("weights.txt", sep="\t")


#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/GO_term_reduction.r --args available_orgs_path=  genomes_path=~/media/jens@margarita/For_CopraRNA2.0/Genomes/ 

# get absolute path  
#path="/home/jens/CopraRNA-git/coprarna_aux/"
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("GO_term_reduction.r","",path)
print(path)


# parameters for term reduction
thres<-0.05		# threshold for enrichment p-Value
max_num=500		# maximal number of members for a term. Removes very broad general terms
min_num=2		# minimal number of members for a term. Renmoves very narrow terms
min_num2=2 		# minimal number of significant genes in a term.
cutoff<-0.5		# overlap threshold for term condensation based on geometric mean of best enrichment p-Value, information content of term and number of significant genes in term
thres_ov<-0.5	# overlap threshold for the number of shared genes between two terms. Proportion of the term with less members.


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

# preset path to required files, path can also be specified as argument
copref_path<-paste(available_orgs_path,"CopraRNA_available_organisms.txt",sep="")
copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	


thres<-as.numeric(thres)
max_num<-as.numeric(max_num)
min_num<-as.numeric(min_num)
min_num2<-as.numeric(min_num2)
cutoff<-as.numeric(cutoff)
thres_ov<-as.numeric(thres_ov)

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

##########
ooi<-readLines("ncrna.fa")
ooi<-gsub(">ncRNA_","",ooi[1])
gbk<-paste(genomes_path,ooi,".gb.gz",sep="")
gbk1<-readLines(gbk)
org<-grep("ORGANISM",gbk1)
org<-gbk1[org]
org<-gsub(".*ORGANISM. ","",org)
org<-gsub(" .*","",org)
gb<-readGenBank(gbk)

cd<-cds(gb)
cd<-GenomicRanges::mcols(cd)

## Now prepare some data.frames
fSym <- cd[,c(3,2)]
fSym<-cbind(seq(1:nrow(cd)),seq(1:nrow(cd)),fSym)
colnames(fSym) <- c("GID","ENTREZID","SYMBOL","GENENAME")

fChr <- cbind(1:nrow(cd), rep(1,nrow(cd)))
fChr<-data.frame(GID=as.integer(fChr[,1]),CHROMOSOME=fChr[,2])
colnames(fChr) <- c("GID","CHROMOSOME")


anno<-enrich_list[[ooi]][[4]]

fGO<-c()
for(i in 1:nrow(fSym)){
	pos<-na.omit(match(fSym[i,3],anno[,13]))
	if(length(pos)>0){
		go<-anno[pos,2]
		go<-strsplit(go, split="; ")[[1]]
		tmp<-cbind(rep(i,length(go)),go)
		fGO<-rbind(fGO,tmp)
	}
}
fGO<-data.frame(GID=as.integer(fGO[,1]),GO=fGO[,2], EVIDENCE=rep("NO",nrow(fGO)))



## Then call the function
makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="blabla",
               genus="A",
               species="B",
               goTable="go")
install.packages("./org.AB.eg.db", repos=NULL)
require(org.AB.eg.db)





#############

term_list<-list()
go<-c("BP","MF","CC")
out_en<-list()
node<-setdiff(tree[[1]][,1],tree[[1]][,2])  # root node
for(ii in 1:3){


ecoGO <- godata("org.AB.eg.db", ont=go[ii])


IC<-c(ecoGO@IC)

	best<-c()
	pval<-c()
	term2<-c()
	inner<-unique(tree[[1]][,1])
	sel<-which(inner==node)
	leaf<-leafs(inner[sel], tree)
	leaf<-tree[[3]][leaf]
	selection<-leaf
	selction<-ord
	term<-c()
	
	for(i in 1:min(length(selection),length(enrich_list))){
		tmp3<-enrich_list[[i]][[1]]
		#tmp<-enrich_list[[i]][[1]][[1]]
		#for(j in 1:3){
		j<-ii
		tmp2<-tmp3[[j]]
		#th<-which(tmp[,"classicFS"]<=thres )
		th2<-which(as.numeric(tmp2[,"classicFS"])<=thres & tmp2[,"Annotated"] > min_num & tmp2[,"Significant"] > min_num2 &  tmp2[,"Annotated"] < max_num)
		term2<-c(term2,tmp2[th2,1])
		pval<-c(pval,tmp2[th2,"classicFS"])
		tmp<-tmp2[th2,2]
		#tmp<-tmp[1:10]
		best<-c(best,tmp2[1:1,2])
		term<-c(term,tmp2[th2,2])
		#}
	}
	k<-which(pval=="< 1e-30")
	if(length(k)>0){
		pval[k]<-"1e-30"
	}
	pval<-as.numeric(pval)
	term<-term[order(pval)]
	term2<-term2[order(pval)]
	pval<-sort(pval)
	count<-table(term)
	term_dup<-which(duplicated(term2))
	numbers<-table(term2)
	term<-term[-term_dup]
	term2<-term2[-term_dup]
	pval<-pval[-term_dup]
	

	best<-unique(best)





term2<-unique(term2)
term_list[[ii]]<-term2
ma<-match(term2,names(numbers))

m2<-match(term2,names(IC))
IC_temp=IC[m2]

na<-which(is.na(IC_temp))

names(IC_temp)[na]<-"GO:1905039"
IC_temp[na]<-0
sim <- mgoSim(term2, term2,
                  semData=ecoGO,
                  measure="Wang",#Jiang
                  combine=NULL)
res<-data.frame(ID=term2,term,pVAL=pval,count=(numbers)[ma], IC=IC_temp)
res[which(is.finite(res$IC)==F),"IC"]<-max(is.finite(res$IC))

full_anno<-c()
for(i in 1:length(enrich_list)){
	tmp3<-enrich_list[[i]][[1]][[ii]]
	if(is.null(tmp3)==F){
		full_anno<-rbind(full_anno,as.matrix(tmp3))
	}
}


gene_cont<-vector("list", length(term2))
names(gene_cont)<-term2
for(i in 1:length(term2)){
	t1<-full_anno[grep(term2[i],full_anno[,1]),"genes"]
	t1<-unlist(strsplit(t1,","))
	gene_cont[[i]]<-t1
}

term_over<-matrix(,length(term2),length(term2))
for(i in 1:(length(term2))){
#print(i)
	for(j in (i):length(term2)){
		t1<-gene_cont[[i]]
		t2<-gene_cont[[j]]
		tmp<-length(intersect(t1,t2))/min(length(t1),length(t2))
		term_over[i,j]<-term_over[j,i]<-tmp
	}
}

number_sig_genes<-rep(NA, length(term2))
names(number_sig_genes)<-term2
for(i in term2){
	tmp<-which(full_anno[,1]==i & full_anno[,"classicFS"] <= thres)
	tmp<-full_anno[tmp,"genes"]
	tmp<-length(unlist(strsplit(tmp, ",")))
	number_sig_genes[i]<-tmp
}


com<-(((-log(as.numeric(res$pVAL))/max(-log(as.numeric(res$pVAL))))**1)*((res$count.Freq/max(res$count.Freq)))*res$IC/max(res$IC[is.finite(res$IC)]))**(1/3)


res$both<-com
res$number_sig_genes<-number_sig_genes

colnames(term_over)<-rownames(term_over)<-term2
ID<-term2	
go1 <- go2 <- overlap <- NULL

ov.df <- as.data.frame(term_over)
ov.df$go1 <- row.names(ov.df)
ov.df <- gather(ov.df, go2, overlap, -go1)
ov.df <- ov.df[!is.na(ov.df$overlap),]				  

	  
				  
by="both"				  


				  
go1 <- go2 <- similarity <- NULL

sim.df <- as.data.frame(sim)
sim.df$go1 <- row.names(sim.df)
sim.df <- gather(sim.df, go2, similarity, -go1)
sim.df <- sim.df[!is.na(sim.df$similarity),]

pos<-na.omit(match(paste(sim.df[,1],sim.df[,2]),paste(ov.df[,1],ov.df[,2])))

sim.df$overlap<-ov.df$overlap[pos]

## feature 'by' is attached to 'go1'
sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")



sim.df$go2 <- as.character(sim.df$go2)


out_en[[ii]]<-list(list(res,sim.df),best)

}

for(i in 1:3){
	tmp<-out_en[[i]][[1]][[1]]
	write.table(tmp, file=paste(go[i],"_enrichment_initial_set.txt",sep=""),sep="\t", quote=F, row.names=FALSE)
}

out_res<-list()

for(j in 1:3){

sim.df<-out_en[[j]][[1]][[2]]
res<-out_en[[j]][[1]][[1]]
best<-out_en[[j]][[2]]
ID<-res$ID
select_fun<-max	


sim.df<-sim.df[order(sim.df$both, decreasing=T),]

GO_to_remove <- character()


    for (i in seq_along(ID)) {
        ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff & sim.df$overlap > thres_ov )

        if (length(ii) < 2)
            next

        sim_subset <- sim.df[ii,]

        jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))

        ## sim.df <- sim.df[-ii[-jj]]
        GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }				  
				  
out2<-res[!res$ID %in% GO_to_remove, ]		
best<-setdiff(best,out2[,2])
if(length(best)>0){
	best<-na.omit(match(best,res[,2]))
	if(length(best)>0){
		out2<-rbind(out2,res[best,])
	}
}
	out2<-out2[order(out2[,"both"], decreasing=T),]	  
	out_res[[go[j]]]<-out2	  
}

names(out_res)<-go


for(ii in go){


out2<-out_res[[ii]]
term<-as.character(out2[,2])
term2<-as.character(out2[,1])
pval<-as.numeric(out2[,"pVAL"])


out<-matrix(,0,4)
colnames(out)<-c("term","org","num","pvalue")
	for(i in 1:length(term)){
		for(j in selection){
			if(is.element(j,names(enrich_list))){
				if(is.null(enrich_list[[j]])==F){
					tmp<-rbind(as.matrix(enrich_list[[j]][[1]][[1]]),as.matrix(enrich_list[[j]][[1]][[2]]),as.matrix(enrich_list[[j]][[1]][[3]]))
					pos<-na.omit(match(term[i], tmp[,2]))
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


ggsave(paste(ii, "_enrichment.pdf",sep=""), height = 3+0.1125*length(tree$tip.label), width = 3.5+0.12*nrow(out2),useDingbats=FALSE)
 

write.table(out_res[[ii]], file=paste(ii, "_enrichment.txt", sep=""), row.names=F,quote=F, sep="\t")



}

	
	

	
	# out<-data.frame(id=as.character(out[,2]),term=as.character(out[,1]), dot_num=as.numeric(out[,3]), pvalue=-1*log10(as.numeric(out[,4])))
	

	# p<-as.numeric(out$pvalue)
	# breaks<-c(0,1.13,2,3,4,5,6)
	# p2<-p
	# for(i in breaks){
		# tmp<-which(p>=i)
		# if(length(tmp)>0){
			# p2[tmp]<-i
		
		# }
	
	# }
	# out$colo<-p2
	# term_map<-out
	# #save(term_map, file="term_map.Rdata")
	# #
	
	# #load("term_map.Rdata")


# col <- c("lightgrey","#7FFFD4","#117893","#34558b","#fd823e","#d13b40" )
# col<-c("#117893","#d13b40" )
# colfunc <- colorRampPalette(col)
# col<-colfunc(6)
# col<-c("grey65",col)
# load("16S_outgroup_tree.Rdata")
# require(forcats)

# #pdf("test2.pdf")

# p <- ggplot(out, aes((term),(id))) + theme_linedraw()
# p + #geom_dotplot(binaxis = "y", stackdir = "center", binpositions="all", binwidth=1, dotsize=0) + 
 # geom_point(aes(size=dot_num, colour=factor(colo))) +
  # theme(legend.text.align = 0,axis.text.x = element_text(angle = 90, size=10, vjust=0.2, hjust=1), axis.text.y = element_text(size=6)) +labs(x="", y="", size="term size", colour="p-Value") + 
  # scale_y_discrete(limits =rev(ord2)) + scale_x_discrete(limits =(term))  + scale_colour_manual(values = col,labels=c("",expression(""<=5%*%10^{-2}),expression(""<=10^{-2}),expression(""<=10^{-3}),expression(""<=10^{-4}),expression(""<=10^{-5}),expression(""<=10^{-6})))
  # # scale_colour_gradientn(colours = colfunc(10), limits=c(0,5),breaks=c(0,1,1.13,2,3,4,5),labels=c(1,expression(10^{-1}),expression(5*10^{-2}),expression(10^{-2}),expression(10^{-3}),expression(10^{-4}),expression(10^{-5})))
 # #dev.off() 
 

# ggsave("test_12.pdf", height = 11, width = 6.5,useDingbats=FALSE)
 