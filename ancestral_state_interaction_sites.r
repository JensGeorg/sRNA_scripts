require(corHMM)
require(seqinr)
require(phangorn)
#require(ggplot2)
#require(ggtree)
#require(ggimage)
#library(NameNeedle)



#call:
#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/ancestral_state_interaction_sites.r 
#R --slave -f /home/jens/media/jens@margarita/CopraRNA-git/coprarna_aux/copraRNA2_find_conserved_sites.r
#path="/home/jens/media/jens@margarita/CopraRNA-git/coprarna_aux/"
#path="/home/jens/CopraRNA-git/coprarna_aux/"


# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("ancestral_state_interaction_sites.r","",path)
#print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	


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
# num<-as.numeric(num)

# sel<-strsplit(sel,split=",")[[1]]
# sel<-as.numeric(sel)


load("16S_rooted_tree.Rdata")
weight<-read.csv("weights.txt", sep="\t")

load("gene_out.Rdata")	

int_sites<-gene_out[[1]]
peak_list<-gene_out[[2]]
col_list<-gene_out[[3]]

# # gain/loss rates of targets based on CopraRNA target predictions

# require(markophylo)
# load("copra_results_all.Rdata")
# load("16S_outgroup_tree.Rdata")

# evo_analysis<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv",sep=",", header=T) 
# evo_analysis<-as.matrix(evo_analysis)


# l<-subtrees(tree)

# tree2<-l[[68]]  #20 eco 16 sal

# m_out<-list()
# for(jj in 1:length(l)){
# tree2<-l[[jj]] 
# all_tars<-c()
# for(i in 1:length(tree2$tip.label)){
	# tmp<-copra_results[[tree2$tip.label[i]]]
	# if(length(tmp)>0){
		# tmp<-as.numeric(gsub(" ","",tmp[,"initial_sorting"]))[1:200]
		# all_tars<-c(all_tars,tmp)
	# }
# }

# all_tars<-unique(all_tars)


# all_targets<-matrix(1,length(all_tars),length(tree2$tip.label))
# colnames(all_targets)<-tree2$tip.label
# rownames(all_targets)<-all_tars
# for(i in 1:length(tree2$tip.label)){
	# tmp<-copra_results[[tree2$tip.label[i]]]
	# if(length(tmp)>0){
		# tmp<-as.numeric(gsub(" ","",tmp[,"initial_sorting"]))[1:200]
		# all_targets[as.character(tmp),tree2$tip.label[i]]<-2
	# }
# }


# fil<-which(rowSums(all_targets)<(ncol(all_targets)+20))


# model1_f<-tryCatch({
    # estimaterates(usertree = tree2, userphyl = all_targets ,
                          # alphabet = c(1, 2), rootprob = "maxlik", 
                          # modelmat = "ARD")
# }, error = function(e) {
    # NULL
# })
# # model1_f <- estimaterates(usertree = tree2, userphyl = all_targets ,
                          # # alphabet = c(1, 2), rootprob = "maxlik", 
                          # # modelmat = "ARD")
						  
# m_out[[jj]]<-model1_f		
			  
# ############################
# }


#selection<-sel

# transform list of sites to matrix
to_table1<-function(x, pos=1, norgs){
	x<-gsub("\\|.*","",x[pos,])
	if(length(x)==0){
		x<-rep(NA,norgs)
	}
	x
}

to_table2<-function(x, pos=1, norgs){
	y<-lapply(x,to_table1, pos=pos, norgs=norgs)
	out<-c()
	for(i in 1:length(y)){
		out<-rbind(out,y[[i]])
	}
	out
}


# tranform the conserved peaks into integers and assign colors for heatmap drawing 
to_number<-function(int_opt, int_sub,homologs, ooi, peaks){
	colo2<-c("#caff70","#7aa9df","#cf97bb" ,"#fede7e","#fc9187","#d3648f","#4c9998","#988377","#576e81","#9fe1e5")
	colo<-c("#799943","#5680b0","#9c6488","#e4c054","#e26a5f","#b43768","#006362","#614736","#294258","#77b1b5")
	ooi_pos<-match(ooi,colnames(int_sub))
	int<-matrix(,nrow(int_opt),ncol(int_opt))
	rownames(int)<-rownames(int_opt)
	colnames(int)<-colnames(int_opt)
	int[]<--1
	int_o<-int
	int_s<-int
	color_vect<-c(-1,0)
	names(color_vect)<-c("#FFFFFF", "#b0b0b0")
	for(i in 1:nrow(int_opt)){
		peaks1<-peaks[[i]]
		exist<-which(homologs[i,]!="")
		if(length(exist)>0){
			int_o[i,exist]<-0
		}
		temp1<-na.omit(unique(c(int_opt[i,])))
		temp2<-na.omit(unique(c(int_sub[i,])))
		temp<-unique(sort(c(temp1,temp2)))
		
		ooi_peak<-int_opt[i,ooi_pos]
		
		if(is.na(ooi_peak)==F){
			ooi_p2<-grep(ooi_peak,peaks1)
			peaks1<-c(peaks1[ooi_p2],peaks1[-ooi_p2])
			
		}
		le<-length(temp)
		if(length(peaks1)>0){
			colo3<-c()
			for(j in 1:length(peaks1)){
				fc <- colorRampPalette(c(colo[j], colo2[j]))
				tmp<-fc(length(peaks1[[j]]))
				names(tmp)<-peaks1[[j]]
				colo3<-c(colo3,tmp)
			}
			no_ex<-which(is.element(colo3, names(color_vect))==F)
			if(length(no_ex)>0){
				tmp<-(max(color_vect)+1):(max(color_vect)+length(no_ex))
				names(tmp)<-colo3[no_ex]
				color_vect<-c(color_vect,tmp)
			}
			colo4<-colo[1:length(peaks1)]
			# for(j in 1:le){
				# tmp<-match(temp[j],names(colo3))
				# tmp<-match(colo3[tmp],names(color_vect))
				# tmp<-color_vect[tmp]
				# pos<-which(int_opt[i,]==temp[j])
				
				# int_o[i,pos]<-tmp
				# pos<-which(int_sub[i,]==temp[j])
				# int_s[i,pos]<-tmp
			# }
			for(j in 1:le){
				tmp<-grep(temp[j],(peaks1))
				#tmp<-match(colo3[tmp],names(color_vect))
				#tmp<-color_vect[tmp]
				pos<-grep(temp[j],int_opt[i,])
				
				int_o[i,pos]<-tmp
				pos<-grep(temp[j],int_sub[i,])
				int_s[i,pos]<-tmp
			}
		}
	}
	out<-list(int_o, int_s, c(c("#FFFFFF", "#b0b0b0"),colo))
	out
}


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




# reference for row positions
evo_analysis<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv",sep=",", header=T) 
evo_analysis<-as.matrix(evo_analysis)


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


fastutr<-read.fasta("utr_seqs.fa")
fastnames<-tolower(names(fastutr)) 


e<-grep("Annotation", colnames(evo_analysis))-1
#genes<-tolower(gsub("\\(.*","",evo_analysis[selection[i],3:e]))
# empty<-which(genes=="")
na<-colnames(evo_analysis)[3:e]
# if(length(empty)>0){
	# genes<-genes[-empty]
	# na<-na[-empty]
# }

out_utrs<-read.fasta("outgroup_UTRs.fa")

# temp2<-na.omit(match(tolower(genes), fastnames))
# write.fasta(fastutr[temp2],file="temp.fa" ,names=na)	
ini<-c()
name<-c()
orgs<-c()
fasta_list<-list()
for(i in 1:length(out_utrs)){
	ids<-names(out_utrs)[i]
	ids<-strsplit(ids,",")[[1]]
	ini<-c(ini,ids[1])
	orgs<-c(orgs,ids[2])
	name<-c(name,ids[3])
	initial_sorting<-as.numeric(ids[1])
	genes<-tolower(gsub("\\(.*","",evo_analysis[initial_sorting,3:e]))
	empty<-which(genes=="")
	na<-colnames(evo_analysis)[3:e]
	if(length(empty)>0){
		genes<-genes[-empty]
		na<-na[-empty]
	}
	temp2<-na.omit(match(tolower(genes), fastnames))
	fast<-fastutr[temp2]
	names(fast)<-na
	fast2<-out_utrs[i]
	names(fast2)<-ids[2]
	fast<-c(fast,fast2)
	fasta_list[[i]]<-fast
	
}





#gene1<-c("galK","rpsQ","sucC","polA","gltA","prpC","purL")
#gene2<-c(1,2,3,4,5,12,13)

anc_res<-vector("list", length(fasta_list))
#names(anc_res)<-gene1
#gggg<-c(15,18,21,23,24)
for(iii in 1:length(fasta_list)){
#for(iii in gggg){


write.fasta(fasta_list[[iii]],file="temp.fa" ,names=names(fasta_list[[iii]]))	

mafft(filename="temp.fa", mode="accurate")
tempf<-read.fasta("ncrna_aligned.fa")
g<-lapply(tempf,paste,collapse="")
g<-unlist(g)
g<-which(duplicated(g))
if(length(g)>0){
	tempf<-tempf[-g]
}
write.fasta(tempf, file.out="ncrna_aligned.fa", names=names(tempf), nbchar=100000)
dat<-read.phyDat("ncrna_aligned.fa", format="fasta", type="DNA")
dm <- dist.ml(dat, model="F81")
treeNJ <- NJ(dm)
fitStart = pml(treeNJ, dat, k=4)
fitJC = optim.pml(fitStart, model="GTR", optGamma=T, rearrangement="stochastic",ratchet.par = list(iter = 5L, maxit = 20L, prop = 1/3),control = pml.control(epsilon = 1e-08, maxit = 10,
trace = 1L))

#tree2<-midpoint(fitJC$tree)
tree2<-root(fitJC$tree,outgroup=names(tempf)[length(tempf)], resolve.root = TRUE)  #psoudomonas as outgroup NC_002516  #galK outgroup "NZ_CP007593"


selection<-as.numeric(ini)

ee<-grep("Annotation", colnames(evo_analysis))-1
homologs<-evo_analysis[selection,3:ee]
org_names<-colnames(evo_analysis)[3:ee]
all_orgs<-colnames(evo_analysis)[3:ee]

int_opt<-to_table2(int_sites[selection], pos=1,norgs=length(org_names))
int_sub<-to_table2(int_sites[selection], pos=2,norgs=length(org_names))

peaks<-peak_list[selection]


ooi<-"NC_000913"
int<-to_number(int_opt,int_sub,homologs,ooi, peaks)

traits<-tree2$tip.label

# variant 2
# trait_list<-list()
# for(i in 1:nrow(int[[1]])){
	
	
	# tmp<-int[[1]][i,]+1
	# tmp2<-int[[2]][i,]+1
	# tmp_out<-rep(0,length(tmp))
	
	# for(j in 0:max(tmp,tmp2)){
		# pos<-which(tmp==j)
		# k<-j
		# if(j>1){
			# pos<-which(tmp>=j)
			# #k<-paste(1,j,sep="&")
		# }
		# if(length(pos)>0){
			# tmp_out[pos]<-k
		# }
	# }		
	# trait_list[[i]]<-data.frame(names(tmp),tmp_out)
# }

# variant 3
trait_list<-list()
for(i in 1:nrow(int[[1]])){
	
	
	tmp<-int[[1]][i,]+1
	tmp2<-int[[2]][i,]+1
	tmp_out<-rep(0,length(tmp))
	
	
	for(j in 1:length(tmp)){
		tmp_out[j]<-tmp[j]
		if(tmp[j]<2 & tmp2[j]>=2){
			tmp_out[j]<-tmp2[j]
		}
		if(tmp[j]<2 & tmp2[j]<2){
			tmp_out[j]<-1
		}
		if(tmp[j]>=2 & tmp2[j]>=2){
			#tmp_out[j]<-paste(sort(c(tmp[j],tmp2[j])),collapse="&")
		}
	}
	pos<-na.omit(match(tree2$tip.label, colnames(int[[1]])))
	#trait_list[[i]]<-data.frame(c(names(tmp)[pos]),c(tmp_out[pos]))
	trait_list[[i]]<-data.frame(c(names(tmp)[pos],setdiff(tree2$tip.label, colnames(int[[1]]))),c(tmp_out[pos],1))
}


recon3<- rayDISC(tree2,data.frame(trait_list[[iii]]),model="ARD", node.states="marginal", root.p="yang",charnum=1)
#recon3<- rayDISC(tree2,data.frame(trait_list[[iii]]),model="ARD", node.states="joint")
#plotRECON(tree2,recon3$states)


orgs<-tree2$tip.label

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

tree3<-tree2
tree3$tip.label<-nam2
anc_res[[iii]]<-list()
anc_res[[iii]][[1]]<-recon3
anc_res[[iii]][[2]]<-tree2
anc_res[[iii]][[3]]<-tree3
anc_res[[iii]][[4]]<-tree
anc_res[[iii]][[5]]<-trait_list[[iii]]
anc_res[[iii]][[6]]<-int[[1]][iii,]+1
anc_res[[iii]][[7]]<-int[[2]][iii,]+1
# pdf(paste("ancestral_states_", gene1[iii],".pdf", sep=""))
	# plotRECON(tree3,recon3$states, piecolors=c(int[[3]][-c(1)]))
# dev.off()


}

names(anc_res)<-name
save(anc_res, file="ancestral_site_recons.Rdata")



dist2root<-function(node, tree){
	root<-setdiff(tree[[1]][,1],tree[[1]][,2])
	count<-0
	while(node!=root){
		count<-count+1
		id1<-which(tree[[1]][,2]==node)	
		node<-tree[[1]][id1,1]		
	}
	return(count)
}



load("ancestral_site_recons.Rdata")
load("copra_results_all.Rdata")
load("target_overlap.Rdata")
load("16S_rooted_tree.Rdata")

cop<-copra_results[[1]]
thres<-0.1
thres2<-0.1

root<-setdiff(tree[[1]][,1],tree[[1]][,2])
b<-tree[[1]][which(tree[[1]][,1]==root),2]
a1<-leafs(b[1],tree)
a2<-leafs(b[2],tree)
node<-b[1]
if(length(a2)>length(a1)){
	node<-b[2]
}
compare_node<-node
overlaps<-a[[1]][[as.character(node)]][[2]]
overlaps<-sort(overlaps, decreasing=T)
overlaps_comp<-a[[1]][[as.character(compare_node)]][[2]]
overlaps_comp<-sort(overlaps_comp, decreasing=T)

overlaps2<-names(overlaps)[which(overlaps>=thres)]
overlaps2_comp<-names(overlaps_comp)[which(overlaps_comp<thres2)]

unique_tars<-intersect(overlaps2,overlaps2_comp)

most_conserved_targets<-cbind(cop[match(overlaps2,cop[,"initial_sorting"]),4], overlaps[which(overlaps>=thres)],overlaps2)



res_table<-matrix(,length(anc_res),8)
colnames(res_table)<-c("name","target_conservation","number_of_sites","sites","organisms_per_site","site_with_highest_lik._at_root","likelihood_at_root","outgroup")





out_list<-list()

for(i in 1:length(anc_res)){
	res_table[i,1]<-names(anc_res)[i]
	pos<-grep(names(anc_res)[i],most_conserved_targets[,1])
	if(length(pos)>0){
		res_table[i,2]<-most_conserved_targets[pos[1],2]
	}

	
	tree<-anc_res[[i]][[3]]
	inner<-unique(tree[[1]][,1])
	root<-setdiff(tree[[1]][,1],tree[[1]][,2])
	states<-anc_res[[i]][[1]]$states
	if(is.null(states)==F){
	node<-sort(inner)
	site_prob<-max(states[which(inner==root),])
	site<-colnames(states)[which(states[which(inner==root),]==max(states[which(inner==root),]))]
	perc<-max(states[which(inner==root),])
	sites<-as.numeric(colnames(states))-1
	nu<-which(sites==0)
	if(length(nu)>0){
		sites<-sites[-nu]
	}
	
	selection<-as.numeric(ini)
	peaks<-peak_list[selection[i]]
	if(as.numeric(site)-1>0){
		site_root<-peaks[[1]][[as.numeric(site)-1]][1]
	} else {
		site_root<-"not_conserved"
	}
	
	sites_id<-c()
	for(j in sites){
		sites_id<-c(sites_id,peaks[[1]][[j]][1])
	}
	names(sites_id)<-sites
	lowest_node<-vector("list",length(sites))
	numb<-c()
	for(j in 1:length(sites)){
		numb[j]<-length(which(anc_res[[i]][[5]][,2]==sites[j]+1))
		tmp<-node[which(states[,as.character(sites[j]+1)]>=0.5)]
		if(length(tmp)>0){
			dists<-unlist(lapply(tmp, dist2root,tree=tree))
			lowest_node[[j]]$lowest_node<-tmp[which(dists==min(dists))]
			lowest_node[[j]]$distance2root<-min(dists)
			orgs_all<-c()
			lowest_node[[j]]$orgs_belonging_to_lowest_node<-list()
			for(jj in 1:length(tmp[which(dists==min(dists))])){				
				orgs<-leafs(tmp[which(dists==min(dists))][j],tree)
				orgs<-tree$tip.label[orgs]
				lowest_node[[j]][[3]][[jj]]<-orgs
				orgs_all<-c(orgs_all,orgs)
			}
			orgs2<-gsub("_.*","",orgs_all)
			lowest_node[[j]]$most_frequent_prefixes<-table(orgs2)
		}
	}
	names(lowest_node)<-sites_id
	out_list[[i]]<-lowest_node
	res_table[i,3]<-length(sites)
	res_table[i,4]<-paste(sites_id,collapse="|")
	res_table[i,5]<-paste(numb,collapse="|")
	if(site_root!="not_conserved"){
		site_root<-paste(site_root,length(which(anc_res[[i]][[5]][,2]==as.numeric(site))),sep="|")
	}
	res_table[i,6]<-site_root
	res_table[i,7]<-perc
	res_table[i,8]<-strsplit(names(out_utrs)[i],",")[[1]][2]
	} else {
		print(i)
	}
}

names(out_list)<-names(anc_res)

res_table<-res_table[order(as.numeric(res_table[,2]),decreasing=T),]


write.table(res_table,file="ancestral_site_reconstruction.txt",sep="\t", row.names=F)
