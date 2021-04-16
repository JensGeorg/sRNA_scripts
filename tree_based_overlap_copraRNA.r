

#path="~/media/jens@margarita/CopraRNA-git/coprarna_aux/"
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggtree))
suppressPackageStartupMessages(require(ggimage))
suppressPackageStartupMessages(require(NameNeedle))
#suppressPackageStartupMessages(require(phangorn))
#suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(ape))

#call:
# R --slave -f ~/media/jens@margarita/CopraRNA-git/coprarna_aux/tree_based_overlap_copraRNA.r --args num=200 outgroup=NC_002516 

num<-200							# number of top targets to be inclued in the analysis
outgroup<-"NC_002516"				# outgroup organism Refseq ID. Outgroup 16S sequence needs to be in input 16S fasta file. In case of "NA" a midpoint rooting is done.
#16S_sequenses<-"16s_sequences.fa" 	# Default for midpoint rooting without outgroup
parts<-2							# Defines the bins for target overlap calculation. 2 equals two bins, i.e <50% and >= 50%.

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("tree_based_overlap_copraRNA.r","",path)


# preset path to required files, path can also be specified as argument
cop_path<-paste(path,"CopraRNA_available_organisms.txt",sep="") # fixxxxxxxxxxxx



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
num<-as.numeric(num)
parts<-as.numeric(parts)

# call mafft for MSA
mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa", mode="fast"){
	if(mode=="accurate"){
		command<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	if(mode=="fast"){
		command<-paste("mafft --retree 2 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	if(mode=="very_fast"){
		command<-paste("mafft --retree 1 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	system(command)
} 

# mafft(filename=16S_sequenses, mode="accurate")
# tempf<-read.fasta("ncrna_aligned.fa")
# write.fasta(tempf, file.out="ncrna_aligned.fa", names=names(tempf), nbchar=100000)
# dat<-read.phyDat("ncrna_aligned.fa", format="fasta", type="DNA")
# dm <- dist.ml(dat, model="F81")
# treeNJ <- NJ(dm)
# fitStart = pml(treeNJ, dat, k=4)
# fitJC = optim.pml(fitStart, model="GTR", optGamma=T, rearrangement="stochastic",ratchet.par = list(iter = 5L, maxit = 20L, prop = 1/3),control = pml.control(epsilon = 1e-08, maxit = 10,trace = 1L))
# fit2_outgroup<-(fitJC$tree)

# if(is.na(outgroup)==F){
	# tree<-root(fitJC$tree, outgroup=outgroup, resolve.root = TRUE)
# } else {
	# tree<-midpoint(fitJC$tree)
# }
# save(tree, file="16S_rooted_tree.Rdata")

load("16S_rooted_tree.Rdata")
load("copra_results_all.Rdata")


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


get_tar_over<-function(copra_results=copra_results, org_list, num=200, weighting=TRUE, int_thres=0.35){
	weight<-read.csv("weights.txt", sep="\t")
	weight<-round(weight,digits=3)*10**3
	tmp<-copra_results[org_list]
	tmp2<-vector("list",length(tmp))
	names(tmp2)<-org_list
	for(i in 1:length(tmp)){		
		tmp3<-tmp[[org_list[i]]]
		if(is.null(tmp3)==F){
		int<-tmp3
		not<-which(int[,org_list[i]]!="")
		int<-int[not,org_list[i]]
		int<-strsplit(int, "\\|")
		int<-as.numeric(do.call(rbind,int)[,3])
		ab<-which(int>int_thres)
		if(length(ab)>0){
			tmp3<-tmp3[-not[ab],]
		}
		tmp3<-tmp3[1:num,"initial_sorting"]
		
		if(weighting==T){
			tmp3<-rep(tmp3,weight[org_list[i],1])
		}
		tmp2[[i]]<-tmp3
		}
	}
	res<-table(unlist(tmp2))
	if(weighting==T){
		res<-res/sum(weight[org_list,])
	} else {
		res<-res/length(org_list)
	}
	res
}


cons_distribution<-function(res, parts=5){
	steps<-1/parts
	out<-vector("list", parts)
	counts<-length(which(res<steps))
	out[[1]]<-names(res)[which(res<steps)]
	for(i in 1:(parts-1)){
		tmp<-length(which(res>=steps*i & res < steps*(i+1)))
		if(i==parts-1){
			tmp<-length(which(res>=steps*i & res <= steps*(i+1)))
		}
		out[[i+1]]<-names(res)[which(res>=steps*i & res < steps*(i+1))]
		if(i==parts-1){
			out[[i+1]]<-names(res)[which(res>=steps*i & res <= steps*(i+1))]
		}
		counts<-c(counts,tmp)
	}
	n1<-c(0,1:(parts-1)*steps)
	n2<-c(steps, 2:parts*steps)
	n<-paste(n1,n2,sep="-")
	names(counts)<-n

	return(list(counts,out))
}


overlap_at_each_node<-function(tree,copra_results=copra_results,parts=5,weighting=TRUE,int_thres=0.35, num=200){
	inner<-unique(tree[[1]][,1])
	over_list<-vector("list", length(inner))
	out<-matrix(,length(inner),parts)
	names(over_list)<-inner
	for(i in 1:length(inner)){
		leaf<-leafs(inner[i], tree)
		leaf<-tree[[3]][leaf]
		over<-get_tar_over(copra_results=copra_results, org_list=leaf, num=num, weighting=weighting, int_thres=int_thres)
		over2<-cons_distribution(over, parts=parts)
		out[i,]<-over2[[1]]
		over_list[[i]]<-list(over2,over)		
	}
	out<-as.data.frame(out)
	out$node<-inner
	return(list(over_list,out))
}


pairwise_identity<-function(tree,weighting=TRUE){
	weight<-read.csv("weights.txt", sep="\t")
	seqs1<-readLines("ncrna.fa")
	na<-grep(">", seqs1)
	seqs<-list()
	for(i in 1:(length(na)-1)){
		seqs[[i]]<-paste(seqs1[(na[i]+1):(na[i+1]-1)],collapse="")
		names(seqs)[i]<-gsub(">ncRNA_","",seqs1[na[i]])
	}
	seqs[[length(na)]]<-paste(seqs1[(na[length(na)]+1):(length(seqs1))],collapse="")
	names(seqs)[length(na)]<-gsub(">ncRNA_","",seqs1[na[length(na)]])
	defaultNeedleParams <- list(MATCH = 5,MISMATCH = -4,GAP = -4,GAPCHAR = "-")
	inner<-unique(tree[[1]][,1])
	topo<-tree[[1]] # tree topology
	branch<-tree[[2]]	
	root<-setdiff(topo[,1],topo[,2])
	inner<-setdiff(inner, root)
	out<-c()
	id_list<-list()
	for(i in 1:length(inner)){
		leaf<-leafs(inner[i], tree)
		leaf<-tree[[3]][leaf]
		ids<-c()
		wei<-c()
		for(j in 1:(length(leaf)-1)){
			seq1<-seqs[[leaf[j]]]
			for(jj in (j+1):length(leaf)){
				seq2<-seqs[[leaf[jj]]]
				compa<-paste(sort(c(leaf[j],leaf[jj])),collapse="_")
				if(is.null(id_list[[compa]])){
					pid <- needles(seq1, seq2, params = defaultNeedleParams)
					pid1 <- strsplit(pid[[2]], "")[[1]]
					pid2 <- strsplit(pid[[3]], "")[[1]]
					id<-length(which(pid1==pid2))/length(pid1)
					id_list[compa]<-id
				} else {
					id<-id_list[[compa]]
				}
				wei<-c(wei,mean(weight[leaf[jj],],weight[leaf[j],]))
				
				ids<-c(ids,id)
			}
			
		}
		if(weighting==T){
			id<-weighted.mean(ids, wei)		
		} else {
			id<-mean(ids)
		}
		out<-c(out, id)
	}      
	out<-data.frame(out, node=inner)
	out
}

# calculate average weighted pairwise sRNA sequence identities at each inner node
out<-pairwise_identity(tree)
save(out, file="pairwise_identity.Rdata")

# bit score

bit<-function(tree){
	weight<-read.csv("weights.txt", sep="\t")
	#seqs<-read.fasta("ncrna.fa",as.string = T)
	seqs1<-readLines("ncrna.fa")
	na<-grep(">", seqs1)
	seqs<-list()
	for(i in 1:(length(na)-1)){
		seqs[[i]]<-paste(seqs1[(na[i]+1):(na[i+1]-1)],collapse="")
		names(seqs)[i]<-gsub(">ncRNA_","",seqs1[na[i]])
	}
	seqs[[length(na)]]<-paste(seqs1[(na[length(na)]+1):(length(seqs1))],collapse="")
	names(seqs)[length(na)]<-gsub(">ncRNA_","",seqs1[na[length(na)]])
	inner<-unique(tree[[1]][,1])
	topo<-tree[[1]] # tree topology
	branch<-tree[[2]]	
	root<-setdiff(topo[,1],topo[,2])
	inner<-setdiff(inner, root)
	out<-c()
	for(i in 1:length(inner)){
		leaf<-leafs(inner[i], tree)
		leaf<-tree[[3]][leaf]
		seq2<-seqs[leaf]
		temp_fasta<-tempfile()
		wtemp<-weight[leaf,]
		wtemp<-wtemp/sum(wtemp)
		names(wtemp)<-leaf
		seq3<-c()
		for(j in 1:length(seq2)){
			seq3<-c(seq3,paste(">",names(seq2)[j],sep=""),seq2[[j]])
		}
		temp_fas<-tempfile()
		writeLines(seq3,con=temp_fas)
		align<-system(paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", temp_fas,seq="") ,intern=T)
		header<-grep(">",align)
		seq4<-vector("list", length(header))
		names(seq4)<-gsub(">","",align[header])
		for(j in 1:(length(header)-1)){
			tmp<-align[(header[j]+1):(header[j+1]-1)]
			tmp<-paste(tmp,collapse="")
			seq4[[j]]<-strsplit(tmp,"")[[1]]
		}
		tmp<-align[(header[length(header)]+1):length(align)]
		tmp<-paste(tmp,collapse="")
		seq4[[length(seq4)]]<-strsplit(tmp,"")[[1]]
		seq4<-do.call(rbind,seq4)
		seq4<-toupper(seq4)
		wtemp<-wtemp*length(leaf)
		new_num<-sum(wtemp)
		score<-c()
		m_score<-c()
		mat<-matrix(,ncol(seq4),4)
		colnames(mat)<-c("A","T","C","G")
		rownames(mat)<-1:ncol(seq4)
		for(j in 1:ncol(seq4)){
			gap<-sum(wtemp[which(seq4[,j]=="-")])
			new_num2<-new_num-gap
			mat[j,1]<-sum(wtemp[which(seq4[,j]=="A")])#/new_num2
			mat[j,2]<-sum(wtemp[which(seq4[,j]=="T")])##/new_num2
			mat[j,3]<-sum(wtemp[which(seq4[,j]=="C")])#/new_num2
			mat[j,4]<-sum(wtemp[which(seq4[,j]=="G")])#/new_num2
			A1<-sum(wtemp[which(seq4[,j]=="A")])/new_num2
			T1<-sum(wtemp[which(seq4[,j]=="T")])/new_num2
			C1<-sum(wtemp[which(seq4[,j]=="C")])/new_num2
			G1<-sum(wtemp[which(seq4[,j]=="G")])/new_num2
			
			en <- min(2,1/(log(2)) * ((4-1)/(2*new_num2)))
			
			tmp_score<-max(log2(4)+sum(min(0,A1*log2(A1),na.rm=T),min(0,T1*log2(T1),na.rm=T),min(0,C1*log2(C1),na.rm=T),min(0,G1*log2(G1),na.rm=T))-en,0)
			max_score<-log2(4)-(min(2,1/(log(2)) * ((4-1)/(2*new_num))))
			score<-c(score,tmp_score)
			m_score<-c(m_score,max_score)
		}
		
		mat<-capture.output(write.table(mat, file=stdout(), sep=" ", quote=F,row.names = T ))
		mat[1]<-paste("P0 ", mat[1],sep="")
		mat<-c("ID any_old_name_for_motif_1","BF species_name_for_motif_1",mat,"XX")
		writeLines(mat, con="MAt.pssm")
		out[i]<-sum(score)/sum(m_score)
	}
	bits<-data.frame(out, node=inner)
	bits
}

bitscore<-bit(tree)

# calculate (weighted) target overlaps at each node
a<-overlap_at_each_node(tree,copra_results=copra_results,parts=parts,weighting=TRUE,num=num)
save(a, file="target_overlap.Rdata")

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
orgs<-tree$tip.label

nam<-c()
	for(i in 1:length(orgs)){
		tnam<-grep(gsub("\\..*","",orgs[i]),copref[,1])
		nam<-c(nam,as.character(copref[tnam[1],2]))
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
fc <- colorRampPalette(c( "#d13b40","#4ec5a5"))
fc <- colorRampPalette(c( "#ffaf12","#8A2BE2"))
fc <- colorRampPalette(c("#117893","#d13b40" ))
tmp<-fc(parts)

dat<-a[[2]] # dataframe to draw nodebars or pies

tree2<-tree
tree2$tip.label<-nam2
d<-data.frame(node=out[,2],num=paste(round(out[,1],digits=2)*100,"%",sep=""))

d2<-data.frame(node=out[,2],score=paste(round(bitscore[,1],digits=2)*100,"%",sep=""))

pdf("traget_overlap_tree.pdf")
p<-ggtree(tree2) %<+% d %<+% d2 + geom_tiplab(size=1.5,hjust=-0.05, align = TRUE, linetype = "dotted", linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6)) + geom_nodelab(aes(label=score), geom = "text", col=1,nudge_x=-0.02, size=2,nudge_y = 0.5) #+ geom_nodelab(aes(label=num), geom = "text", col=1,nudge_x=-0.02, size=2,nudge_y = 0.5)
bars <- nodebar(dat, cols=1:parts, color=tmp,position='dodge',alpha=1)
p1<-inset(p, bars, width=0.051, height=1)
p1
dev.off()

## plot node reference tree 
pdf("tree_node_reference.pdf")
plot(tree2, cex=0.3)
nodelabels(cex=0.3, frame="none")
dev.off()