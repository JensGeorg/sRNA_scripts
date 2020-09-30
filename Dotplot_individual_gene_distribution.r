
require(ggtree)
require(ggplot2)

#call:
#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/Dotplot_individual_gene_distribution.r --args sel=866,1494,1472,1466,4617,1465,3449,3752,2057,2490,2967,1075,3358,14360,6799,4208,5037,3046,4682,14540 
#R --slave -f /home/jens/media/jens@margarita/CopraRNA-git/coprarna_aux/copraRNA2_find_conserved_sites.r
#path="/home/jens/media/jens@margarita/CopraRNA-git/coprarna_aux/"
#path="/home/jens/CopraRNA-git/coprarna_aux/"


# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("Dotplot_individual_gene_distribution.r","",path)
#print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	

num=200 # maximal CopraRNA rank considered to be a positive prediction
ooi<-"NC_000913" # organisms used for annotations


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


sel<-strsplit(sel,split=",")[[1]]
sel<-as.numeric(sel)
#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/conserved_site_on_subset.r --args ooi=NC_000911 sel=10836,12682,12166,10810,11659,13432,12099,12247,12267,12453,10453,13292,12882,11396,12443,2171,11690,10726,11535,1203,12702,1515

ref<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv",sep=",")

load("16S_rooted_tree.Rdata")
load("gene_out.Rdata")
load("target_overlap.Rdata")	
load("copra_results_all.Rdata")

int_sites<-gene_out[[1]]
peak_list<-gene_out[[2]]
col_list<-gene_out[[3]]



leafs<-function(node, tree){  # finds leafs belonging to an inner node
	topo<-tree[[1]]
	branch<-tree[[2]]
	names(node)<-1
	out<-node
	edge<-c()
	nn<-node
	temp<-c()
	
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


orgs<-tree$tip.label

nam<-c()
	for(i in 1:length(orgs)){
		tnam<-grep(gsub("\\..*","",orgs[i]),copref[,1])
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

nam2<-paste(nam2,orgs,sep="_")

cop<-copra_results[[ooi]]
cop2<-copra_results[[ooi]]
node<-setdiff(tree[[1]][,1],tree[[1]][,2])  # root node
if(sum(a[[1]][[as.character(node)]][[1]][[1]])==0){
	pos<-which(tree[[1]][,1]==node)
	pos<-tree[[1]][pos,2]
	no<-c()
	no[1]<-length(leafs(pos[1],tree))
	no[2]<-length(leafs(pos[2],tree))
	node<-pos[which(no==max(no))[1]]
}
overlaps<-a[[1]][[as.character(node)]][[2]]
overlaps<-sort(overlaps, decreasing=T)
overlaps2<-names(overlaps)


#hh<-cbind(cop2[match(overlaps2,cop2[,"initial_sorting"]),ooi], overlaps,overlaps2)
e<-which(colnames(ref)=="Annotation")-1
pos<-match(sel,as.numeric(overlaps2))
hh<-cbind(ref[sel,3:e], overlaps[pos],as.numeric(overlaps2[pos]))
nam<-c()
for(i in 1:nrow(hh)){
	tmp<-gsub(".*\\(","",hh[i,1:(ncol(hh)-3)])
	tmp<-gsub("\\|.*","",tmp)
	na<-which(tmp!="N/A" & tmp!= "")
	if(length(na)==0){
		tmp<-gsub("\\(.*","",hh[i,1:(ncol(hh)-3)])
		na<-which(tmp!="N/A" & tmp!= "")
		if(length(na)==0){
			tmp<-hh[i,ncol(hh)]
			na<-1
		}
	} 
	
	nam<-c(nam,tmp[na[1]])
}	

inner<-unique(tree[[1]][,1])
sel2<-which(inner==node)
leaf<-leafs(inner[sel2], tree)
leaf<-tree[[3]][leaf]
selection<-leaf



out_num<-matrix(,nrow(hh),length(selection))
colnames(out_num)<-selection
rownames(out_num)<-nam




# for(i in 1:nrow(hh)){
	# if(is.na(nam[i])){
		# tmp<-NA
		# j<-1
		# while(is.na(tmp) & j<=length(selection)){
			# tmp<-match(hh[i,3], copra_results[[selection[j]]][,"initial_sorting"])
			# j<-j+1
		# }
		# if(is.na(tmp)==F){
			# tmp<-copra_results[[selection[j-1]]][tmp,]
			# e<-grep("Annotation", names(tmp))
			# tmp<-tmp[4:(e-1)]
			# tmp<-tmp[which(tmp!="")]
			# gene<-gsub(".*\\(","",tmp)
			# gene<-unlist(lapply(gene, function(x){return(strsplit(x,"\\|")[[1]][1])}))
			# locus<-gsub("\\(.*","",tmp)
			# ge<-which(gene!="N/A")
			# if(length(ge)>0){
				# nam[i]<-gene[ge[1]]
			# } else {
				# nam[i]<-locus[1]
			# }
		# }
	# }
# }

na<-which(is.na(nam))
if(length(na)>0){
	nam<-nam[-na]
	hh<-hh[-na,]
}

org<-c()
genes<-c()
ranks<-c()
color<-c()

for(i in 1:nrow(hh)){
	for(j in 1:length(selection)){
		
		or<-grep(selection[j],nam2)
		org<-c(org,nam2[or])
		
		
		tmp<-match(hh[i,ncol(hh)], as.numeric(copra_results[[selection[j]]][,"initial_sorting"]))
		if(is.na(tmp)==F){
		print(j)
		}
		ranks<-c(ranks,tmp)
		genes<-c(genes,nam[i])
		
		int<-int_sites[[as.numeric(gsub(" ","",hh[i,ncol(hh)]))]]
		peaks_all<-peak_list[[as.numeric(gsub(" ","",hh[i,ncol(hh)]))]]
		colo4<-col_list[[as.numeric(gsub(" ","",hh[i,ncol(hh)]))]]
		
		if(length(peaks_all)>0){
		
		int_tmp<-match(selection[j], colnames(int))
		int_tmp<-int[,int_tmp]
		if(is.na(int_tmp[1])==F){
			int1<-int_tmp[1]
			int1<-strsplit(int1,"\\|")[[1]][1]
		} else {
			int1<-int_tmp[2]
			int1<-strsplit(int1,"\\|")[[1]][1]
		}
		if(is.na(int1)){
			color<-c(color, "#b0b0b0")
		} else {
			tmp_col<-match(int1, names(colo4))
			if(is.na(tmp_col)){
				color<-c(color, "#b0b0b0")
			} else {
				color<-c(color, colo4[tmp_col])
			}
		}
		} else {
			color<-c(color, "#b0b0b0")
		}
		
		
		
	}
}
ranks2<-ranks
ranks2[which(ranks<=num)]<-num-as.numeric(as.character(ranks2[which(ranks<=num)]))
ranks2[which(ranks>num)]<-0.01

mi<-which(is.na(ranks2)==F & is.na(color)==T) 
if(length(mi)>0){
	color[mi]<-1
}
cols<-unique(color)
names(cols)<-1:length(cols)
for(i in 1:length(cols)){
	color[grep(cols[i],color)]<-i
}

out3<-data.frame(org=as.character(org),genes=genes,ranks=ranks2, color=color, ori_ranks=ranks)


d = fortify(tree)
d = subset(d, isTip)
ord<-with(d, label[order(y, decreasing=T)])
ord2<-c()
for(i in ord){
	p<-grep(i, nam2)
	ord2<-c(ord2,nam2[p])
}

p <- ggplot(out3, aes((genes),(org))) + theme_linedraw()
p + 
geom_point(aes(size=ranks, colour=factor(color))) + scale_y_discrete(limits =rev(ord2)) + scale_x_discrete(limits =(nam)) + scale_colour_manual(values = cols) + labs(x="", y="", size="prediction rank") +
 theme(axis.text.x = element_text(angle = 90, size=10, vjust=0.5, hjust=1, face="italic"))
 
 
 ggsave("dotplot_conservation_selected_genes.pdf",height = 3+0.1125*length(tree$tip.label), width = 3.5+0.25*nrow(hh),useDingbats=FALSE)


# cd CopraRNA-git/coprarna_aux/	
 # +
  # theme(axis.text.x = element_text(angle = 45, size=8, vjust=1, hjust=1), axis.text.y = element_text(size=6)) +labs(x="", y="", size="prediction rank") + 
  # scale_y_discrete(limits =rev(ord2)) + scale_x_discrete(limits =(term))  + scale_colour_manual(values = col)
 