# Wrapper script for a CopraRNA target prediction / evolutionary analysis
# Starting point is a GLASSgo sRNA homolog prediction fasta file

# extract synteny network, visualize syntenies sRNA distribution and conserved sNA and promter elements
## Requirements: Syntney (https://github.com/CyanolabFreiburg/Syntney)

### create network
python3 ~/media/jens@margarita/Syntney/Syntney.py -i ~/media/jens@margarita/Copra2_paper/Glassgo/RyhB/RyhB_ref2.fa  -o ~/media/jens@margarita/Copra2_paper/Glassgo/RyhB/ -n cys -r off -d ~/Syntney_db/synt.db -c ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r  -s ~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py
python3 ~/media/jens@margarita/Syntney/Syntney.py -i [path to GLASSgo fasta]  -o [path to output folder] -n cys -r off -d [path to Syntney database] -c [./Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r]  -s [./Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py]

### annotate and reduce network
R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/network_modifications.r --args working_directory="" only_sig_nodes=TRUE max_synt=20 thres_anno=0.025 thres_edge=0.025

### select representative homologs, draw sRNA/Synteny distribution tree, draw promoter/sRNA weblogo, draw linear syntenies of representative candidates
R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/sRNA_promoter_conservation.r --args filename=~/media/jens@margarita/Syntney/testfiles/inputForJens.fasta  synteny_window=5000 script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney/new.db



# Script runs after a standard CopraRNA run





# Make predictions for all organisms (if not already done by the standard call)
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/join_pvals_coprarna_2.r --args ooi_only=FALSE prediction_on_subset=TRUE	outlier_removal=TRUE maxorgs=15")
R --slave -f ~/CopraRNA-git/coprarna_aux/join_pvals_coprarna_2.r --args ooi_only=FALSE prediction_on_subset=TRUE outlier_removal=TRUE maxorgs=15


# calculate target overlap and sRNA sequence identities at each internal node 
# tree is calculated based on the 16S sequences of the organisms 
# an outgroup needs to be added to the 16S fasta file or the tree is midpoint rooted

suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(seqinr))

## call mafft for MSA
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

## make reference tree, outgroup or midpoint rooted
### add 16S sequence of an appropriate outgroup organisms to the copraRNA generated 16S fasta file (Fasta header should be Refseq ID)
ribo_sequenses="16s_sequences_outgroup.fa"
outgroup<-"NC_002516"	
mafft(filename=ribo_sequenses, mode="accurate")
tempf<-read.fasta("ncrna_aligned.fa")
write.fasta(tempf, file.out="ncrna_aligned.fa", names=names(tempf), nbchar=100000)
dat<-read.phyDat("ncrna_aligned.fa", format="fasta", type="DNA")
dm <- dist.ml(dat, model="F81")
treeNJ <- NJ(dm)
fitStart = pml(treeNJ, dat, k=4)
fitJC = optim.pml(fitStart, model="GTR", optGamma=T, rearrangement="stochastic",ratchet.par = list(iter = 5L, maxit = 20L, prop = 1/3),control = pml.control(epsilon = 1e-08, maxit = 10,trace = 1L))
fit2_outgroup<-(fitJC$tree)

if(is.na(outgroup)==F){
	tree<-root(fitJC$tree, outgroup=outgroup, resolve.root = TRUE)
} else {
	tree<-midpoint(fitJC$tree)
}
save(tree, file="16S_rooted_tree.Rdata")

# ## plot tree 
# pdf("tree_node_reference.pdf")
# plot(tree, cex=0.5)
# nodelabels(cex=0.5, frame="none")
# dev.off()

## calculate target overlap and average pairwise sRNA identities at each internal node of the rference tree
## output: "traget_overlap_tree.pdf" + "tree_node_reference.pdf" to identify the the IDs of the internal nodes
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/tree_based_overlap_copraRNA.r --args num=200 outgroup=NC_002516 16S_sequenses=16s_sequences_outgroup.fa")
R --slave -f ~/media/jens@margarita/CopraRNA-git/coprarna_aux/tree_based_overlap_copraRNA.r --args num=200 outgroup=NC_002516

## find conserved targets at root or another selectable node. Nodes are displayed in the tree_node_reference.pdf
## find targets that are unique to a given node in ccomaprison with anoter node

load("copra_results_all.Rdata")
load("target_overlap.Rdata")
load("16S_rooted_tree.Rdata")

cop<-copra_results[[1]]
thres<-0.8
thres2<-0.15
node<-156 #156 #161
compare_node<-161
overlaps<-a[[1]][[as.character(node)]][[2]]
overlaps<-sort(overlaps, decreasing=T)
overlaps_comp<-a[[1]][[as.character(compare_node)]][[2]]
overlaps_comp<-sort(overlaps_comp, decreasing=T)

overlaps2<-names(overlaps)[which(overlaps>=thres)]
overlaps2_comp<-names(overlaps_comp)[which(overlaps_comp<thres2)]

unique_tars<-intersect(overlaps2,overlaps2_comp)

most_conserved_targets<-cbind(cop[match(overlaps2,cop[,"initial_sorting"]),4], overlaps[which(overlaps>=thres)],overlaps2)
most_conserved_targets

unique_tars<-cbind(cop[match(unique_tars,cop[,"initial_sorting"]),4],unique_tars)
unique_tars


node<-101 #156 #161
compare_node<-117
thres<-0.1
thres2<-0.4
overlaps<-a[[1]][[as.character(node)]][[2]]
overlaps<-sort(overlaps, decreasing=T)
overlaps_comp<-a[[1]][[as.character(compare_node)]][[2]]
overlaps_comp<-sort(overlaps_comp, decreasing=T)
overlaps2<-names(overlaps)[which(overlaps<thres)]
overlaps2_comp<-names(overlaps_comp)[which(overlaps_comp>=thres2)]
lost<-setdiff(overlaps2_comp,names(overlaps))
lost<-cbind(cop[match(lost,cop[,"initial_sorting"]),4],lost)
lost

## percentage of target conservation at a selected node

node<-156
node_cons<-a[[2]][which(a[[2]][,3]==node),]
more_or_equal_50<-node_cons[2]/sum(node_cons[1:2])
less_tan_50<-node_cons[1]/sum(node_cons[1:2])
out<-c(more_or_equal_50,paste(node_cons[2],sum(node_cons[1:2]),sep="/"),less_tan_50,paste(node_cons[1],sum(node_cons[1:2]),sep="/"))
names(out)<-c("more_or_equal_50%","more_or_equal_50%","less_than_50%","less_than_50%")
out



# do functional enrichment
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/functional_enrichment.r --args num=200")

# apply automatic GO-term complexity reduction
# output PDF-files and tables can be used for final manual selection
# parameters for term reduction
thres=0.05		# threshold for enrichment p-Value
max_num=500		# maximal number of members for a term. Removes very broad general terms
min_num=2		# minimal number of members for a term. Renmoves very narrow terms
min_num2=2 		# minimal number of significant genes in a term.
cutoff=0.5		# overlap threshold for term condensation based on geometric mean of best enrichment p-Value, information content of term and number of significant genes in term
thres_ov=0.5	# overlap threshold for the number of shared genes between two terms. Proportion of the term with less members.

system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/GO_term_reduction.r --args thres=0.05 max_num=500 min_num=2 min_num2=2 cutoff=0.5 thres_ov=0.5")


# Details of the enrichment
load("enrichment.Rdata")
#number of orgs with annotation
length(enrich_list)-length(which(unlist(lapply(enrich_list, is.null))))

#members of a given term for organisms belonging to a a given internal node


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


node_enrich_summary<-function(tree, enrich_list,node=95, term="chemotaxis",num=20, max_num=500, min_num=5){
	inner<-unique(tree[[1]][,1])
	sel<-which(inner==node)
	leaf<-leafs(inner[sel], tree)
	leaf<-tree[[3]][leaf]
	selection<-leaf
	out<-matrix(,0,6)
	colnames(out)<-c("locus_tag","gene","enrich_pval","organism","CopraRank","description")
	
	anno<-c()
	for(i in 1:length(selection)){
		tmp1<-enrich_list[[selection[i]]][[3]]
		if(is.null(tmp1)==F){
			anno<-rbind(anno,as.matrix(tmp1))
		}
	}	
	enri<-c()
	for(i in 1:length(selection)){
		tmp<-enrich_list[[selection[i]]][[1]][[1]]
		tmp2<-enrich_list[[selection[i]]][[1]][[2]]
		tmp3<-enrich_list[[selection[i]]][[1]][[3]]
		enri<-rbind(enri,tmp,tmp2,tmp3)
	
	}
	
			tmp<-grep(term,enri[,1])
			if(length(tmp)>0){
				#tmp<-tmp[which(tmp[,"Annotated "]] < max_num & tmp[,"Annotated "]] >= min_num & tmp[,"Significant "]] > 0),]
				#tmp<-tmp[1:num,]
				#pos<-na.omit(match(term,tmp[,2]))
				#if(length(pos)!=0){			
					genes<-strsplit(enri[tmp,"genes"],split=",")
					le<-unlist(lapply(genes,length))
					genes<-unlist(genes)
					pval<-rep(enri[tmp,"classicFS"],le)
					ranks<-unlist(strsplit(enri[tmp,"ranks"],split=","))
					#org<-selection[i]
					#org<-grep(org,nam2)
					#org<-rep(nam2[org],length(genes))
					pos<-match(toupper(genes),toupper(anno[,ncol(anno)-1]))
					gene<-as.character(anno[pos,9])
					des<-as.character(anno[pos,12])
					org<-as.character(anno[pos,5])
					tmp1<-cbind(genes,gene,pval,org,ranks,des)
					out<-rbind(out,tmp1)
				
				#}			
			}
		out
	}
genes<-node_enrich_summary(tree, enrich_list, node=161, term="GO:0006099")
sort(table(node_enrich_summary(tree, enrich_list, node=156, term="GO:0051536")[,2]), decreasing=T)

# Next manually review enrichment files and merge them into a final file (just delete un-wanted rows). Do final visualization of enrichment
## produces "final_enrichment.pdf"
enrich<-read.csv("selected_enrichment.txt", sep="\t")
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/plot_enrichment.r --args enrich_file=selected_enrichment.txt")




# Do site conservation analysis on subset of selected targets and generate 3D interaction probability landscape
## select your favorite targets for visualization
## define the selection by the "initial_sorting"-ID from the organism specific result files or the  

###Spot42
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/conserved_site_on_subset.r --args sel=866,1494,1472,1466,4617,1465,3449,3752,2057,2490,2967,1075,3358,14360,6799,4208,5037,3046,4682,14540")
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/conserved_site_on_subset.r --args sel=2490,2967,3358,14360,6799,4208,5037,3046,4682,14540")
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/conserved_site_on_subset.r --args sel=2611,3366,4160,3126,4114,4002,965,4525,3524")
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/conserved_site_on_subset.r --args sel=3404,2950,4683")
###RyhB


###sel<-c(1184,685,3501,1185,609,2050,2006,4147,4373,3663,1347,1430,2548,2623,2067,3017,4617,2620,3642,1230,12205,1426,626,10915,5170,1932,4104) # RyhB

# Do conservation dotplots for selected targets
## Output is the "dotplot_conservation_selected_genes.pdf" file
system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/Dotplot_individual_gene_distribution.r --args sel=866,1494,1472,1466,4617,1465,3449,3752,2057,2490,2967,1075,3358,14360,6799,4208,5037,3046,4682,14540,3407,4064,4757,3350,3669,3143,2950,4321,5752")

# Do ancestral site reconstruction on subset of selected targets and plot the respective trees
## requires a fasta file with outgroup UTRs for each selected UTR (outgroup_UTRs.fa), the fasta header should be the initial_sorting"-ID , the outgroup refseq ID or another organism identifier and an optionally description separated by an ",".
## the outgroup_UTRs.fa defines which UTRs are analyzed. The UTRs needed to be part of the above subset analysis.

system("R --slave -f /home/jens/CopraRNA-git/coprarna_aux/ancestral_state_interaction_sites.r")

## plot trees with most likely ancestral states

system("R --slave -f ~/media/jens@margarita/CopraRNA-git/coprarna_aux/plot_ancestral_states.r")




