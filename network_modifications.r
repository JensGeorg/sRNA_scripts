
#CALL:
#R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/network_modifications.r --args working_directory="" only_sig_nodes=TRUE max_synt=20 thres_anno=0.025 thres_edge=0.025

working_directory<-getwd()
only_sig_nodes<-TRUE
max_synt<-20
# thresholds for including nodes or edges into the annotated network
thres_anno<-0.025
thres_edge<-0.025


args <- commandArgs(trailingOnly = TRUE) 

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
}
only_sig_nodes<-as.logical(only_sig_nodes)
max_synt<-as.numeric(max_synt)
thres_anno<-as.numeric(thres_anno)
thres_edge<-as.numeric(thres_edge)

setwd(working_directory)
d<-dir()

network_file<-d[grep("_Network.txt",d)[1]]
anno_file<-d[grep("_Network_Annotation.txt",d)[1]]
cluster_file<-d[grep("_Network_Cluster.txt",d)[1]]
synt_table<-d[grep("_Synteny_Table.txt",d)[1]]


synt<-read.csv(synt_table,sep="\t", header=T)
network<-read.csv(network_file,sep=",", header=T)
anno<-read.csv(anno_file,sep="\t", header=T)

pos<-na.omit(match(network[,1],anno[,1]))
PageRank<-rep(NA, nrow(anno))
anno<-cbind(anno,PageRank)

# Cluster with a score below the threshold are flagged to be drawn transparently
anno[pos,"PageRank"]<-network[,"PageRank"]
low<-which(anno[,"PageRank"]<thres_anno)
filtered_names<-as.character(anno[,2])
transparent<-rep("no",nrow(anno))


if(length(low)>0){
	filtered_names[low]<-""
	transparent[low]<-"yes"
}
anno<-cbind(anno,filtered_names,transparent)

cluster<-read.csv(cluster_file,sep="\t", header=T) # file containing the existing syntenys in "cluster" annotation

upstream<-unlist(lapply(as.character(cluster[,2]), strsplit, split=","))
downstream<-unlist(lapply(as.character(cluster[,3]), strsplit, split=","))
al<-c(upstream,downstream)
al<-al[-grep("sRNA",al)]
al<-sort(table(al),decreasing=T)

# create matrix that counts co-appearance of protein clusters in the input syntenys
# cluster are sorted decreasingly by the number of co-appearances 
clustab<-matrix(0,length(al),length(al))  
colnames(clustab)<-rownames(clustab)<-names(al)

for(i in 1:nrow(cluster)){
	up<-unlist(lapply(as.character(cluster[i,2]), strsplit, split=","))
	down<-unlist(lapply(as.character(cluster[i,3]), strsplit, split=","))
	up<-up[-grep("sRNA",up)]
	down<-down[-grep("sRNA",down)]
	if(length(up)>0 & length(down) >0){
	for(j in 1:length(up)){
		for(jj in 1:length(down)){
			clustab[up[j],down[jj]]<-clustab[up[j],down[jj]]+1
			clustab[down[jj],up[j]]<-clustab[down[jj],up[j]]+1
		}
	}
	if(length(up)>1){
		for(j in 1:(length(up)-1)){
			for(jj in (j+1):length(up)){
				clustab[up[j],up[jj]]<-clustab[up[j],up[jj]]+1
				clustab[up[jj],up[j]]<-clustab[up[jj],up[j]]+1
			}
		}
	}
	if(length(down)>1){
		for(j in 1:(length(down)-1)){
			for(jj in (j+1):length(down)){
				clustab[down[j],down[jj]]<-clustab[down[j],down[jj]]+1
				clustab[down[jj],down[j]]<-clustab[down[jj],down[j]]+1
			}
		}
	}
	}
}

# greedy method to combine co-appearing clusters in synteny families.
# clusters are joined to a query cluster if >=50% of its total apearance counts
# are through the query cluster. The procedure starts with the
# family reference cluster (the custer with the highest appearance in the family)
# all clusters joined to a family serve as new query clusters as long as no new connected members
# are detected. Finally all clusters belonging to this network are removed from the count table, 
# making the procedure a hard clustering method.


clus_list<-list()
i<-1
while(ncol(clustab)>1){
	tmp<-c(colnames(clustab)[1])
	memb<-colnames(clustab)[1]
	inves<-c()
	while(length(memb)>0){
		count<-al[memb[1]]
		inves<-c(inves,memb[1])
		memb2<-which(clustab[,memb[1]]>=1)#>count/2)	
		count1<-clustab[memb2,memb[1]]
		memb2<-rownames(clustab)[memb2]
		count2<-al[memb2]
		count3<-count1/count2
		memb2<-names(which(count3>=0.60))
		memb2<-setdiff(memb2,tmp)
		tmp<-c(tmp,memb2)
		memb<-setdiff(tmp,inves)
	}
	del<-na.omit(match(tmp, colnames(clustab)))
	clustab<-clustab[-del,-del]
	clus_list[[i]]<-tmp
	i<-i+1
}


score<-rep(NA,length(clus_list))
names(score)<-1:length(score)

for(i in 1:length(clus_list)){
	tmp<-match(clus_list[[i]],anno[,1])
	tmp<-anno[tmp,"PageRank"]
	tmp<-sum(as.numeric(tmp))
	score[i]<-tmp
}
score<-sort(score,decreasing=T)

clus_list<-clus_list[as.numeric(names(score))]
names(clus_list)<-1:length(clus_list)


# detecting the organism subsets represented in the synteny families
# annotating the clusters with the included organisms
# create table that matches the unique sRNA IDs with synteny families

orgs_accesion<-rep(NA, nrow(anno))
orgs<-rep(NA, nrow(anno))
anno<-cbind(anno,orgs,orgs_accesion)
org_list<-list()
id_fam<-c()

for(i in 1:length(clus_list)){
	tmp<-clus_list[[i]]
	tmp_orgs<-c()
	tmp_ids<-c()
	for(j in tmp){
		
		on<-grep(paste(j,",",sep=""), cluster[,2])
		tw<-grep(paste(j,",",sep=""), cluster[,3])
		on<-c(on,tw)
		if(length(on)>0){
			on<-cluster[on,1]
			tmp_ids<-c(tmp_ids,on)
			on<-gsub("_.*","",on)
			on2<-org_names<-unlist(lapply(on, function(x){
					return(gsub(" .*","",synt[match(x,synt[,"accesion"])[1],"organism"]))	
				}))
			on2<-names(sort(table(org_names),decreasing=T))
			anno[match(j,anno[,1]),"orgs_accesion"]<-paste(on, collapse=",")
			anno[match(j,anno[,1]),"orgs"]<-paste(on2, collapse=",")
			tmp_orgs<-c(tmp_orgs, on)
		}		
	}
	tmp_orgs<-unique(tmp_orgs)
	tmp_ids<-unique(tmp_ids)
	tmp_ids<-cbind(rep(names(clus_list)[i],length(tmp_ids)),tmp_ids)
	id_fam<-rbind(id_fam,tmp_ids)
	org_names<-unlist(lapply(tmp_orgs, function(x){
		return(gsub(" .*","",synt[match(x,synt[,"accesion"])[1],"organism"]))	
	}))
	org_names<-sort(table(org_names),decreasing=T)
	org_list[[i]]<-org_names
}



cluster_id<-rep(NA,nrow(anno))
for(i in 1:nrow(anno)){
	tmp<-which(is.na(lapply(clus_list, function(x){match(as.character(anno[i,1]),x)}))==F)
	if(length(tmp)>0){
		cluster_id[i]<-tmp
	}

}
anno<-cbind(anno,cluster_id)


##############
id_fam<-c()
l2<-lapply(clus_list,paste0,"_")

mat<-c()
for(i in 1:length(l2)){
	tmp<-unlist(clus_list[i])
	tmp<-cbind(paste0(tmp,"_"), rep(names(clus_list)[i],length(tmp)))
	mat<-rbind(mat,tmp)
}

for(i in 1:nrow(cluster)){
	print(i)
	tmp<-unlist(strsplit(as.character(cluster[i,2:3]),","))
	tmp<-paste0(tmp,"_")
	count<-mat[match(tmp,mat[,1]),2]
	
	
	# count<-lapply(tmp, function(x){
		# return(grep(paste0(x,"_"),l2))
	# })
	# count<-unlist(count)
	if(length(count)>0){
		pr<-anno[match(tmp,paste0(anno[,1],"_")),"PageRank"]
		names(pr)<-count
		pr<-na.omit(pr)
		count<-sort(table(count),decreasing=T)#[1]
		count2<-rep(0,length(count))
		names(count2)<-names(count)
		for(j in 1:length(count)){
			tm<-sum(pr[which(names(pr)==names(count)[j])])
			count2[j]<-tm
		}
		count2<-sort(count2,decreasing=T)
		id_fam<-rbind(id_fam, c(cluster[i,1],names(count2)[1],paste(names(count2),collapse=",")))
	} else{
		id_fam<-rbind(id_fam, c(cluster[i,1],NA,NA))
	}	
}
id_fam<-id_fam[,c(2,1,3)]

# edges with a combined weight below the given threshold are drawn transparently 
combined.weight<-network[,"PageRank"]*network[,"connection.weight"]

transparent<-rep("no",nrow(network))
low<-which(combined.weight<thres_edge)


if(length(low)>0){
	transparent[low]<-"yes"
}

network<-cbind(network, combined.weight, transparent)


# identification if the cluster is predominately 5' (left) or 3' (right) of the reference gene
cluster_id<-rep(NA,nrow(network))
orientation<-rep(NA,nrow(network))
#position_relative_to_sRNA<-rep(NA,nrow(anno))

#anno<-cbind(anno,position_relative_to_sRNA)

for(i in 1:nrow(anno)){
	tmp<-which(as.character(network[,1])==as.character(anno[i,1]))
	left<-length(grep(paste0(as.character(anno[i,1]),","),cluster[,2]))
	right<-length(grep(paste0(as.character(anno[i,1]),","),cluster[,3]))
	anno[i,"position_relative_to_sRNA"]<-"left"
	if(right>left){
		anno[i,"position_relative_to_sRNA"]<-"right"
	}
	if(length(tmp)>0){
		cluster_id[tmp]<-as.character(anno[i,"cluster_id"])
		orientation[tmp]<-as.character(anno[i,"position_relative_to_sRNA"])
	}
}
network<-cbind(network,cluster_id, orientation)

write.table(network, file=paste("post_processed_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno, file=paste("post_processed_",anno_file, sep=""),sep="\t", row.names=F, quote=F)


# drawing the detected families as un connected sub-networks
# if only_sig_nodes is TRUE end-nodes below the given threshold are removed
# from the network.

check_sig<-function(network){
	out<-setdiff(network[,1],network[,2])
	tr<-c()
	pos<-c()
	for(i in out){
		tmp_pos<-match(i, anno[,1])
		if(anno[tmp_pos,"transparent"]=="yes"){
			tr<-c(tr, i)
		}
	}
	for(i in tr){
		pos<-c(pos,match(i,network[,1]))
	}			
	while(length(pos)>0){
		network<-network[-pos,]
		out<-setdiff(network[,1],network[,2])
		tr<-c()
		pos<-c()
		for(i in out){
			tmp_pos<-match(i, anno[,1])
			if(anno[tmp_pos,"transparent"]=="yes"){
				tr<-c(tr, i)
			}
		}
		for(i in tr){
			pos<-c(pos,match(i,network[,1]))
		}	
	}
	out<-setdiff(network[,2],network[,1])
	tr<-c()
	pos<-c()
	out<-setdiff(out,"sRNA")
	for(i in out){
		tmp_pos<-match(i, anno[,1])
		if(anno[tmp_pos,"transparent"]=="yes"){
			tr<-c(tr, i)
		}
	}
	for(i in tr){
		pos<-c(pos,match(i,network[,2]))
	}	
	if(length(pos)>0){
		network<-network[-pos,]
	}
	network
}



not_connected2sRNA<-function(tmp){
	ids<-unique(tmp[,1])
	connected<-c()
	srna<-grep("sRNA",tmp[,2])
	while(length(srna)>0){
		pos<-srna[1]		
		pos_del<-pos
		#connected<-c(connected,tmp[pos,1])
		while(length(pos)>0){
			connected<-c(connected,tmp[pos,1])
			tmp2<-tmp[pos,1]
			pos2<-c()
			tmp<-tmp[-pos_del,]			
			for(j in tmp2){
				pos2<-c(pos2, na.omit(which(tmp[,2]==j)))
			}
			pos_del<-pos2
			pos<-unique(pos2)
			
		}	
		srna<-grep("sRNA",tmp[,2])
	}
	connected<-unique(connected)
	not_connected<-setdiff(ids,connected)
	not_connected
}

surf<-function(network, inp){
	pos<-c()
	con<-c()
	for(i in 1:length(inp)){
		pos2<-which(network[,1]==inp[i])
		tar<-network[pos2,2]
		sr<-match("sRNA", tar)
		if(is.na(sr)==F){
			con<-c(con,inp[i])
		} else {
			names(tar)<-paste(rep(names(inp[i]),length(tar)),rep(inp[i],length(tar)),sep=",")
			pos<-c(pos,tar)			
		}
	}	
	if(length(con)>0){
		con<-paste(names(con),con,sep=",")
		con<-sub(",","",con)
		con<-strsplit(con,split=",")
		return(con)
	}
	return(surf(network, pos))
}

uid<-names(clus_list)

network2<-c()

for(i in 1:min(length(uid),max_synt=max_synt)){
	tmp<-which(network[,"cluster_id"]==uid[i])
	tmp<-network[tmp,]
	
	if(only_sig_nodes==TRUE){
		tmp<-check_sig(tmp)
	}
	
	notcon<-not_connected2sRNA(tmp)
	if(length(notcon)>0){
		notcon<-notcon[which(tmp[match(notcon,tmp[,1]),"PageRank"]>=thres_anno)]
		if(length(notcon)>0){
			notcon<-notcon[order(tmp[match(notcon,tmp[,1]),"PageRank"], decreasing=T)]
			path<-surf(network, notcon)
			if(length(path)>1){
				#find path with highest overlap to sub network
				ov<-unlist(lapply(path, function(x){
					tmp<-length(match(x,c(tmp[,1],tmp[,2])))
				}))
				ov<-which(ov==max(ov))
				path<-path[[ov]]
				if(length(path)>1){
					# length of path to srna
					le<-unlist(lapply(path, length))
					le<-which(le==min(le))
					path<-path[[le]]
					if(length(path)>1){
						path<-path[[1]]
					}
				}
			} else {
				path<-path[[1]]
			}
			path<-c(path,"sRNA")
			for(j in 1:(length(path)-1)){
				te<-which(network[,1]==path[j] & network[,2]==path[j+1])
				tmp<-rbind(tmp,network[te,])
			}
		}
	}
	
	if(nrow(tmp)>0){
		tmp_ids<-unique(c(as.character(tmp[,1]),as.character(tmp[,2])))
		tmp_anno<-anno[na.omit(match(tmp_ids,anno[,1])),]
		tmp_anno[,1]<-paste(tmp_anno[,1],i,sep="_")
		
		tmp[,1]<-paste(tmp[,1],i,sep="_")
		tmp[,2]<-paste(tmp[,2],i,sep="_")
		network2<-rbind(network2,tmp)
		anno<-rbind(anno,tmp_anno)
	}
}



write.table(network2, file=paste("post_processed2_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno, file=paste("post_processed2_",anno_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(id_fam, file="ids2synt_families.txt",sep="\t", row.names=F, quote=F)





tnet<-which(network2[,"transparent"]=="yes")
if(length(tnet)>0){
	network2<-network2[-tnet,]
}

write.table(network2, file=paste("post_processed_wo_trans_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno, file=paste("post_processed_wo_trans_",anno_file, sep=""),sep="\t", row.names=F, quote=F)









