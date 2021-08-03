
#CALL:
#R --slave -f  ~/media/jens@margarita/sRNA_scripts/network_modifications2.r --args script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney_db/synt.db filename=~/media/jens@margarita/Copra2_paper/Glassgo/Spot42/Spot42_rev.fa working_directory=~/media/jens@margarita/Copra2_paper/Glassgo/Spot42/ only_sig_nodes=TRUE thres_anno=0.001 thres_edge=0.001
#R --slave -f  ~/media/jens@margarita/sRNA_scripts/network_modifications2.r --args script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/synteny/synt.db filename=~/synteny/rsaE_glassgo.fa working_directory=~/synteny/ only_sig_nodes=TRUE thres_anno=0.001 thres_edge=0.001



# dependencies: 
require(stringi)
require(phangorn)
require(ggtree)
require(MCL)
#require(RCy3)
#genbank_groper_sqliteDB_ver00.py
#CALL:
options(scipen=999)
#filename<-file('stdin', 'r') # result fasta file from GLASSgo
#filename<-"~/For_CopraRNA2.0/OxyS/IsaR1/clear/IsaR1clear.txt"

filename<-"/home/jens/media/jens@margarita/Staph_enrichment/01_GLASSgo_Results/01_GLASSgo_Results/01_GLASSgo_Results/01_GLASSgo_Results/GLASSgo_output_HG001_01716.fa"
filename<-"/home/jens/media/jens@margarita/Copra2_paper/Glassgo/RyhB/RyhB_ref2.fa"
#filename<-"~/media/jens@margarita/Copra2_paper/Glassgo/Spot42/Spot42_rev.fa"
script_path<-"~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
#script_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
#db_path<-"/media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db"
db_path<-"~/synteny/synt_rRNA.db"
working_directory<-getwd()
only_sig_nodes<-TRUE
max_synt<-1000
# thresholds for including nodes or edges into the annotated network
thres_anno<-0.001
thres_edge<-0.001
up<-150
down<-50
threads<-40
rRNA_existence_threshold<-0.6

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
up<-as.numeric(up)
down<-as.numeric(down)

setwd(working_directory)
d<-dir()
network_file<-d[grep("_Network.txt",d)[1]]
anno_file<-d[grep("_Network_Annotation.txt",d)[1]]
cluster_file<-d[grep("_Network_Cluster.txt",d)[1]]
synt_table<-d[grep("_Synteny_Table.txt",d)[1]]


synt<-read.csv(synt_table,sep="\t", header=T)
network<-read.csv(network_file,sep=",", header=T)
anno<-read.csv(anno_file,sep="\t", header=T)


split_glassgo<-function(x){
	tmp<-strsplit(x[1], ">")[[1]]
	if(length(tmp)>2){
		tmp2<-strsplit(tmp[length(tmp)],"p.c.VAL:")[[1]]
		tmp2<-gsub(";.*","",tmp2[length(tmp2)])
		tmp<-paste(tmp[2], "p.c.VAL:", tmp2[length(tmp2)], sep="")
	}else {
			tmp<-tmp[2]
		}
	tmp<-strsplit(tmp, ":")[[1]]
	id<-tmp[1]
	taxid<-tmp[length(tmp)]
	identity<-strsplit(tmp[length(tmp)-1],"-taxID")[[1]][1]
	tmp2<-strsplit(tmp[2]," ")[[1]]
	tmp3<-paste(tmp2[2:(length(tmp2)-1)],collapse=" ")
	name<-strsplit(tmp3,",")[[1]][1]
	tmp2<-strsplit(tmp2,"-")[[1]]
	a<-tmp2[1]
	b<-as.numeric(tmp2[2])
	a<-strsplit(a,"c")[[1]]
	strand<-"+"
	a<-as.numeric(a[length(a)])
	if(b<a){
		st<-b
		en<-a
		strand<-"-"
	}else{
		st<-a
		en<-b
	}
	ID<-paste(id,a,sep="_")
	out<-c(id, strand, st,en,name,ID)
	out
}

export_ncRNA_coordinates<-function(x){ 
	header_row <- grep(">", x)
	headers <- as.character(x[header_row])
	first_line<-grep(":",headers[1])
	seqs<-as.character(x[header_row+1])
	if(length(first_line)==0){
		headers<-headers[2:length(headers)]
		seqs<-seqs[2:length(seqs)]
	}
	
	tmp<-do.call(rbind,lapply(headers,split_glassgo))
	tmp<-cbind(tmp,headers,seqs)
	colnames(tmp)<-c("Accesion_number", "Strand","start","end","name","ID","Full_header","sequence")
	tmp
}


locus_tag2org<-function(out2){
	tag<-c()
	org<-c()
	for(i in 1:length(out2)){
		temp_tag<-out2[[i]][,5]
		temp_org<-rep(names(out2)[i], length(temp_tag))
		tag<-c(tag,temp_tag)
		org<-c(org,temp_org)
	}
	out<-cbind(tag,org)
	out
}


rand_extension<-function(x, Accession){
	ra<-stri_rand_strings(x,length=4,  pattern = "[A-Za-z0-9]")
	temp<-paste(Accession,ra,sep="_")
	#temp<-gsub("\"","",temp)
	
	temp
}

get_prot_fasta3<-function(out, filen){
  fasta<-c()
  for(i in 1:length(out)){
    for(j in 1:nrow(out[[i]])){
        if(is.na(out[[i]][j,6])==F){
        temp<-as.character(out[[i]][j,6])
		na<-as.character(out[[i]][j,5])
		na<-gsub("\\\"","",na)
         na<-paste(">",na,sep="")
		 temp<-c(na,temp)
         fasta<-c(fasta,temp)
        }
      }
    }
  write.table(fasta, file=filen, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}


# identify homologous proteins using CDhit
cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T, threads=2){
	tempf<-tempfile()
	wd<-getwd()
	di<-paste(wd,"/", "psi_out",sep="")
	dir.create(di)
	if(psi==T){
		inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -d 50 -o ", tempf, " -c ", thres, sep="")
	 }
	if(psi==F){
		inp<-paste("cd-hit -i ", fasta,  " -d 50 -o ",tempf, " -c ",  thres ," -n 2", " -aL 0.6", " -T ", threads, sep="")
	}
	 # print(inp)
	 system(inp, intern=T)
	 cd<-paste(tempf, ".clstr", sep="")
	 cd<-readLines(cd)
	 cd<-as.character(cd)
	 cd<-gsub("\t"," ", cd)
	 cd
}

proc_cdhit<-function(x){ 
  clustlist<-list()
  numb<-grep(">Cluster", x)
  for(i in 1:length(numb)){
    if(i<length(numb)){
      end<-numb[i+1]-1
    }
    if(i==length(numb)){
      end<-length(x)
    }
    temp<-x[(numb[i]+1):end]
    temp<-gsub(".*aa, >","",temp)
    temp<-gsub("\\.\\.\\..*","",temp)
    clustlist[[i]]<-temp
  }
  clustlist	
}



remove_overlapping_homologs<-function(coor, over_thres=0.5){
	orgs<-unique(coor[,1])
	over<-c()
	for(i in orgs){
		tmp<-which(coor[,1]==i)
		if(length(tmp)>1){
			coordina<-function(coor){
				out<-as.numeric(coor[3]):as.numeric(coor[4])
				out
			}
			coordi<-apply(coor[tmp,],1,coordina)
			if(is.matrix(coordi)==T){
				coordi<-as.list(data.frame(coordi))
			}
			names(coordi)<-coor[tmp,"ID"]
			
			j<-1
			while(length(coordi)>1 & j<length(coordi)){
			
				tmp<-c()
				for(jj in (j+1):length(coordi)){
					l1<-length(coordi[[j]])
					l2<-length(coordi[[jj]])
					l<-min(l1,l2)
					int<-length(intersect(coordi[[j]],coordi[[jj]]))/l
					names(int)<-paste(names(coordi)[j],names(coordi)[jj],sep=";")
					tmp<-c(tmp, int)
					
				}
				tmp2<-which(tmp>over_thres)
				
				if(length(tmp2)>0){
					na<-strsplit(names(tmp)[tmp2],";")[[1]][2]
					over<-c(over, na)
					coordi<-coordi[-(tmp2+1)]
					
				}
				j<-j+1
			}			
		}	
	}
	if(length(over)>0){
		over<-na.omit(match(over, coor[,"ID"]))
		if(length(over)>0){
			coor<-coor[-over,]
		}
	}
	coor
}

make_oneline_fasta<-function(x){
	header <- grep(">", x)
	header<- c(header,length(x)+1)
	out<-c()
	for(i in 1:(length(header)-1)){
		tmp<-paste(x[(header[i]+1):(header[i+1]-1)],collapse="")
		out<-c(out, x[header[i]],tmp)
	}
	out
}


# Execute functions
fasta<-readLines(filename)
closeAllConnections()
#fasta<-as.character(fasta)
fasta<-make_oneline_fasta(fasta)
coor<-export_ncRNA_coordinates(fasta)
while(length(which(duplicated(coor[,"ID"])))>0){
	coor<-remove_overlapping_homologs(coor)
}

orgs<-(coor[,1])
command<-paste("python3 ", script_path, " -s ", db_path, " -rs ", paste(orgs, collapse=" "))
print(command)
refseq<-system(command, intern=T)


odd<-seq(1,length(refseq)-1,by=2)
refseq2<-refseq[odd]
refseq2<-strsplit(refseq2,split="is")
refseq2<-do.call(rbind,refseq2)
refseq2<-gsub("refseq ID of  ","",refseq2)
refseq2<-gsub("!","",refseq2)
refseq2<-gsub(" ","",refseq2)
colnames(refseq2)<-c("identifier","fin")

coor<-cbind(coor,refseq2)
save(coor, file="coor.Rdata")

#16S RNA
orgs<-unique(coor[,1])
command<-paste("python3 ", script_path, " -s ", db_path, " -rRNA ", paste(orgs, collapse=" "))
#print(command)
rRNA<-system(command, intern=T)

startp<-grep(">",rRNA)[1]
rRNA<-rRNA[startp:length(rRNA)]

empty<-grep("No 16srRNA found",rRNA)
empty2<-which(rRNA=="")
empty<-c(empty,empty2)
if(length(empty)>0){
	rRNA<-rRNA[-c(empty-1,empty)]
}


# check number of input sequences with coreesponding rRNA
orgs<-grep(">", rRNA)
orgs2<-gsub(">","",rRNA[orgs])


rRNA2<-c()
count<-0
for(i in 1:length(orgs)){
	tmp<-grep(orgs2[i], coor[,1])
	if(length(tmp)>0){
		count<-count+length(tmp)
		rRNA3<-matrix(,length(tmp),2)
		for(j in 1:length(tmp)){
			rRNA3[j,1]<-coor[tmp[j],"Full_header"]
			rRNA3[j,2]<-rRNA[orgs[i]+1]
		}
		rRNA2<-rbind(rRNA2,rRNA3)
	}
}

removed2<-c()
if(count/nrow(coor)>=rRNA_existence_threshold){
	orgs<-grep(">", rRNA)
	orgs<-gsub(">","",rRNA[orgs])
	removed<-c()
	for(i in 1:nrow(coor)){
		tmp<-na.omit(match(coor[i,1],orgs))
		if(length(tmp)==0){
			removed<-c(removed,i)
		}
	}
	if(length(removed)>0){
		removed2<-coor[removed,"ID"]
		coor<-coor[-removed,]
	}	
} 
if(length(removed2)>0){
	removed2<-cbind(removed2, rep("no_16S_rRNA",length(removed2)))
}


# create bed style file for sequence extraction


if(nrow(coor)>1000){
	jobs1<-nrow(coor)%/%1000
	jobs1<-max(2,jobs1)
	#jobs<-jobs*1000
	jobs<-nrow(coor)%/%jobs1
	rest<-nrow(coor)-(jobs*jobs1)
	jobs<-rep(jobs,jobs1)
		if(rest>0){
			jobs[1:rest]<-jobs[1:rest]+1
		}
	count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
	count_vect2<-cumsum(jobs)
} else {
	jobs<-nrow(coor)
	jobs1<-1
	rest<-0
	count_vect1<-1
	count_vect2<-nrow(coor)
}



bed3<-c()
for(j in 1:length(count_vect1)){
	bed<-c()

	for(i in count_vect1[j]:count_vect2[j]){
		if(coor[i,2]=="+"){
			tmp<-c(coor[i,1],as.numeric(coor[i,3])-up,as.numeric(coor[i,4])+down, coor[i,2])
		} else {
			tmp<-c(coor[i,1],as.numeric(coor[i,3])-down,as.numeric(coor[i,4])+up, coor[i,2])
		}
		bed<-rbind(bed,tmp)
	}

	tmpf<-tempfile()
	tmpf<-"temp.txt"

	bed2<-c()
	for(i in 1:nrow(bed)){
		bed2<-paste(bed2,paste(bed[i,],collapse="@"),sep=" ")
	}
	bed2<-gsub("^ ","",bed2)
	#write.table(bed, file=tmpf, quote=F, sep=";", row.names=F, col.names=F)
	#command<-paste("python3 ", script_path, " -s ", db_path, " -pdna ", tmpf, sep="")
	command<-paste("python3 ", script_path, " -s ", db_path, " -pdna ", bed2,sep="")
	bed2<-system(command, intern=T)
	bed3<-c(bed3,bed2)
}
sta<-grep(">", bed3)[1]
bed2<-bed3[sta:length(bed3)]

combined_seqs<-c()

#rRNA2[,1]<-gsub(">","",rRNA2[,1])
#rRNA2[,1]<-gsub(":.*","",rRNA2[,1])

count<-0

for(i in 1:nrow(coor)){
	tmp<-c()
	r<-na.omit(which(paste(">",coor[i,1],sep="")==rRNA)[1])
	if(length(r)>0){
		if( rRNA[r+1]!="No 16sRNA sequence found"){
			tmp<-c(tmp, coor[i,"ID"],coor[i,"Full_header"],coor[i,"sequence"],rRNA[r+1])
			if(coor[i,2]=="+"){
				na<-paste(">",coor[i,1],":",as.numeric(coor[i,3])-up,"-",as.numeric(coor[i,4])+down,sep="")
			}
			if(coor[i,2]=="-"){
				na<-paste(">",coor[i,1],":c",as.numeric(coor[i,4])+up,"-",as.numeric(coor[i,3])-down,sep="")
			}
			pr<-na.omit(which(bed2==na)[1])
			if(length(pr)>0){
				#print(i)
				tmp<-c(tmp,bed2[pr+1])		
				combined_seqs<-rbind(combined_seqs,tmp)
			} else{
				count<-count+1
				print(i)
			}
		} 
	} 
}

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


tree_weights<-function(tree, method="clustal"){
	tip<-Ntip(tree)
	node<-Nnode(tree)
	di<-dist.nodes(tree)
	root<-setdiff(tree[[1]][,1],tree[[1]][,2])
	out<-vector("list",tip)
	names(out)<-1:tip
	for(i in 1:tip){
		temp<-numeric(0)
		t1<-i
		while(t1!=(tip+1)){
			t1<-tree[[1]][which(tree[[1]][,2]==t1),1]
			temp<-c(temp,t1)
		}
		out[[i]]<-temp
	}
	count<-table(unlist(out))
	le<-0
	for(i in 1:length(count)){
		t1<-tree[[1]][which(tree[[1]][,2]==as.numeric(names(count)[i])),1]
		le<-c(le,di[t1,as.numeric(names(count)[i])])
	}
	if(method=="clustal"){
		count<-le/count
	}
	if(method=="copra"){
		names(le)<-names(count)
		count<-le/2
		
	}
	weight<-numeric()
	for(i in 1:tip){
		t1<-di[as.numeric(names(out)[i]),out[[i]][1]]
		po<-match(out[[i]],names(count))
		weight<-c(weight,sum(count[po],t1))
		
	}
	names(weight)<-tree$tip.label
	weight<-weight/sum(weight)
	weight
}


# weights based on 16S rDNA
fasta16<-c()
for(i in 1:nrow(combined_seqs)){
	fasta16<-c(fasta16,paste(">",combined_seqs[i,1], sep=""))
	fasta16<-c(fasta16, combined_seqs[i,4])
}

temp_fas<-tempfile()
temp_fas2<-tempfile()
writeLines(fasta16,con=temp_fas)
align16<-system(paste("mafft --thread ", threads," --retree 2 --maxiterate 0 --quiet --inputorder ", temp_fas, " > ", temp_fas2,seq=""))

ribo<-read.phyDat(temp_fas2, format="fasta", type="DNA")
dm <- dist.ml(ribo, model="F81")
fitJC<- upgma(dm)
weight16<-tree_weights(fitJC, method="clustal")
if(length(unique(weight16))==1){
	if(unique(weight16)=="NaN"){
		weight16<-rep(1,length(weight16))
	}
}
dist_ham_ribo<-dist.hamming(ribo)
dist_ham_ribo<-as.matrix(dist_ham_ribo)


# weights based on sRNA
fasta_sRNA<-c()
for(i in 1:nrow(combined_seqs)){
	fasta_sRNA<-c(fasta_sRNA,paste(">", combined_seqs[i,1], sep=""))
	fasta_sRNA<-c(fasta_sRNA, combined_seqs[i,3])
}

temp_fas<-tempfile()
temp_fas2<-tempfile()
writeLines(fasta_sRNA,con=temp_fas)
align_sRNA<-system(paste("mafft --thread ", threads," --retree 2 --maxiterate 0 --quiet --inputorder ", temp_fas, " > ", temp_fas2,seq=""))
srna<-read.phyDat(temp_fas2, format="fasta", type="DNA")
dm <- dist.ml(srna, model="F81")
fitJC<- upgma(dm)
weight_sRNA<-tree_weights(fitJC, method="clustal")
if(length(unique(weight_sRNA))==1){
	if(unique(weight_sRNA)=="NaN"){
		weight_sRNA<-rep(1,length(weight_sRNA))
		}
}
dist_ham_srna<-dist.hamming(srna)
dist_ham_srna<-as.matrix(dist_ham_srna)

# combined weights
weight_comb<-apply(cbind(weight16,weight_sRNA), 1, function(x){return(exp(mean(log(x))))})
#weight_comb<-weight_comb/max(weight_comb)
weight_comb<-weight_comb/sum(weight_comb)


save(weight_comb, file="phylogenetic_weights.Rdata")
save(combined_seqs, file="combined_seqs.Rdata")
dist_objects<-list(dist_ham_srna,dist_ham_ribo)
save(dist_objects, file="dist_objects.Rdata")

cluster<-read.csv(cluster_file,sep="\t", header=T) # file containing the existing syntenys in "cluster" annotation
#ex<- na.omit(match(cluster[,1],coor[,"ID"]))
#coor<-coor[ex,]
cluster<-cluster[na.omit(match(names(weight_comb),cluster[,1])),]

pos<-na.omit(match(network[,1],anno[,1]))
PageRank<-rep(NA, nrow(anno))
anno<-cbind(anno,PageRank)

# Cluster with a score below the threshold are flagged to be drawn transparently
anno[pos,"PageRank"]<-network[,"Sum.of.branches"]
low<-which(anno[,"PageRank"]<thres_anno)
filtered_names<-as.character(anno[,2])
transparent<-rep("no",nrow(anno))


if(length(low)>0){
	filtered_names[low]<-""
	transparent[low]<-"yes"
}
anno<-cbind(anno,filtered_names,transparent)

upstream<-unlist(lapply(as.character(cluster[,2]), strsplit, split=","))
downstream<-unlist(lapply(as.character(cluster[,3]), strsplit, split=","))
al<-c(upstream,downstream)
al<-al[-grep("sRNA",al)]
clus<-unique(al)
#al<-sort(table(al),decreasing=T)

# create matrix that counts co-appearance of protein clusters in the input syntenys
# cluster are sorted decreasingly by the number of co-appearances 
clustab<-matrix(0,length(clus),length(clus))  
colnames(clustab)<-rownames(clustab)<-clus

for(i in 1:nrow(cluster)){
	up<-unlist(lapply(as.character(cluster[i,2]), strsplit, split=","))
	down<-unlist(lapply(as.character(cluster[i,3]), strsplit, split=","))
	up<-up[-grep("sRNA",up)]
	down<-down[-grep("sRNA",down)]
	id<-cluster[i,1]
	pos<-match(id,names(weight_comb))
	if(is.na(pos)){
		print(i)
	}
	id_weight<-weight_comb[pos]
	if(length(up)>0 & length(down) >0){
		for(j in 1:length(up)){
			for(jj in 1:length(down)){
				clustab[up[j],down[jj]]<-clustab[up[j],down[jj]]+id_weight
				clustab[down[jj],up[j]]<-clustab[down[jj],up[j]]+id_weight
			}
		}
	} 
	if(length(up)>1){
			for(j in 1:(length(up)-1)){
				for(jj in (j+1):length(up)){
					clustab[up[j],up[jj]]<-clustab[up[j],up[jj]]+id_weight
					clustab[up[jj],up[j]]<-clustab[up[jj],up[j]]+id_weight
				}
			}
	}
	if(length(down)>1){
			for(j in 1:(length(down)-1)){
				for(jj in (j+1):length(down)){
					clustab[down[j],down[jj]]<-clustab[down[j],down[jj]]+id_weight
					clustab[down[jj],down[j]]<-clustab[down[jj],down[j]]+id_weight
				}
			}
	}
}
#clustab<-clustab/max(clustab)

#al<-colSums(clustab)
#al<-apply(clustab,1,max)
#al<-sort(al, decreasing=T)

#clustab<-clustab[names(al),names(al)]

######
for(i in 1:ncol(clustab)){
	clustab[,i]<-clustab[,i]/max(clustab[,i])
}

cl<-mcl(clustab, addLoops=F)

if(cl=="An Error occurred at iteration 1"){
	clus_list<-rep(1,nrow(clustab))
} else{
	clus_list<-cl$Cluster

}
out<-list()
for(i in 1:max(clus_list)){
	tmp<-which(clus_list==i)
	out[[i]]<-colnames(clustab)[tmp]
}
emp<-which(unlist(lapply(out, length))==0)
if(length(emp)>0){
	out<-out[-emp]
}
clus_list<-out
###

# greedy method to combine co-appearing clusters in synteny families.
# clusters are joined to a query cluster if >=50% of its total apearance counts
# are through the query cluster. The procedure starts with the
# family reference cluster (the custer with the highest appearance in the family)
# all clusters joined to a family serve as new query clusters as long as no new connected members
# are detected. Finally all clusters belonging to this network are removed from the count table, 
# making the procedure a hard clustering method.


# c2<-clustab
# clus_list<-list()
# i<-1
# while(ncol(clustab)>1){
	# tmp<-c(colnames(clustab)[1])
	# memb<-colnames(clustab)[1]
	# inves<-c()
	# while(length(memb)>0){
		# #count<-al[memb[1]]
		# inves<-c(inves,memb[1])
		# #memb3<-clustab[,memb[1]]/al[colnames(clustab)]
		# #memb2<-which(clustab[,memb[1]]==0.01 | memb3 > 0.1)#>count/2)
		# na<-colnames(clustab)
		# memb2<-which(clustab[,memb[1]]==al[na])
		# #memb2<-which(clustab[,memb[1]]>0.5)
		# # count1<-clustab[memb2,memb[1]]
		# # memb2<-rownames(clustab)[memb2]
		# # count2<-al[memb2]
		# # count3<-count1/count2
		 # memb2<-colnames(clustab)[memb2]
		# memb2<-setdiff(memb2,tmp)
		# tmp<-c(tmp,memb2)
		# memb<-setdiff(tmp,inves)
	# }
	# del<-na.omit(match(tmp, colnames(clustab)))
	# clustab<-clustab[-del,-del]
	# clus_list[[i]]<-tmp
	# i<-i+1
	# if(is.vector(clustab)==T){
		# break
	# }
# }

# # sort synteny family clusters by the sum of the normalized counts of their member protein clusters

# c2<-clustab
# min_con<-7/nrow(cluster)*5
# synt_families<-function(clustab, min_connectivity=0){
	# clus_list<-list()
	# i<-1
	# while(ncol(clustab)>1){
		# tmp<-c(colnames(clustab)[1])
		# memb<-colnames(clustab)[1]
		# inves<-c()
		# while(length(memb)>0){
			# #count<-al[memb[1]]
			# inves<-c(inves,memb[1])
			# #memb3<-clustab[,memb[1]]/al[colnames(clustab)]
			# #memb2<-which(clustab[,memb[1]]==0.01 | memb3 > 0.1)#>count/2)
			# na<-colnames(clustab)
			# memb2<-c()
			# for(j in 1:nrow(clustab)){
				# if(clustab[memb[1],j]==al[memb[1]] | clustab[memb[1],j]==al[rownames(clustab)[j]] ){
					# if(clustab[memb[1],j] >= min_connectivity){
						# memb2<-c(memb2,rownames(clustab)[j])
					# }
				# }
			# }
			# #memb2<-which(clustab[,memb[1]]==al[na])
			# #memb2<-which(clustab[,memb[1]]>0.5)
			# # count1<-clustab[memb2,memb[1]]
			# # memb2<-rownames(clustab)[memb2]
			# # count2<-al[memb2]
			# # count3<-count1/count2
			# #memb2<-colnames(clustab)[memb2]
			# memb2<-setdiff(memb2,tmp)
			# tmp<-c(tmp,memb2)
			# memb<-setdiff(tmp,inves)
		# }
		# del<-na.omit(match(tmp, colnames(clustab)))
		# clustab<-clustab[-del,-del]
		# clus_list[[i]]<-tmp
		# i<-i+1
		# if(is.vector(clustab)==T){
			# break
		# }
	# }
	# clus_list
# }

# clus_list<-synt_families(clustab, min_connectivity=min_con)

# singletons<-which(unlist(lapply(clus_list, length))==1)
# cl<-clus_list[singletons]
# clus_list<-clus_list[-singletons]

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

#membership

memb<-matrix(0,length(clus),length(clus_list))
rownames(memb)<-clus
colnames(memb)<-names(clus_list)
fam_vec<-c()

for(i in 1:nrow(memb)){
	tmp_sum<-sum(clustab[,clus[i]])
	fam_vec<-c(fam_vec,grep(paste0(clus[i],"$"), clus_list))
	for(j in 1:ncol(memb)){
		clus_sum<-sum(clustab[clus_list[[j]],clus[i]])
		memb[i,j]<-clus_sum/tmp_sum
	}
}
memb<-cbind(memb,fam_vec)
write.table(memb, file="synteny_family_membership.txt", sep="\t", quote=F)

# left<-c()
# for(i in 1:length(cl)){
	# tmp<-which(memb[cl[[i]],]==max(memb[cl[[i]],]))[1]
	# if(memb[cl[[i]],tmp]>0){
		# clus_list[[tmp]]<-c(clus_list[[tmp]],cl[[i]])
	# } else{
		# left<-c(left,cl[[i]])
	# }
# }

# if(length(left)>0){
	# clustab<-c2[left,left]
	# al<-apply(clustab,1,max)
	# al<-sort(al, decreasing=T)
	# clustab<-clustab[names(al),names(al)]
	# clus_list<-synt_families(clustab, min_connectivity=0)



# }

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
	tmp_ids<-cbind(rep(names(clus_list)[i],length(tmp_ids)),as.character(cluster[tmp_ids,1]))
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
	#print(i)
	tmp<-unlist(strsplit(c(as.character(cluster[i,3]),as.character(cluster[i,2])),","))
	tmp<-paste0(tmp,"_")
	count<-na.omit(mat[match(tmp,mat[,1]),2])
	
	
	# count<-lapply(tmp, function(x){
		# return(grep(paste0(x,"_"),l2))
	# })
	# count<-unlist(count)
	if(length(count)>0){
		pr<-anno[match(as.character(tmp),paste0(anno[,1],"_")),"PageRank"]
		names(pr)<-count
		pr<-as.numeric(na.omit(pr))
		
		count<-sort(table(count),decreasing=T)#[1]
		count2<-rep(0,length(count))
		names(count2)<-names(count)
		for(j in 1:length(count)){
			tm<-sum(pr[which(names(pr)==names(count)[j])])
			count2[j]<-tm
		}
		count2<-sort(count2,decreasing=T)
		id_fam<-rbind(id_fam, c(as.character(cluster[i,1]),names(count2)[1],paste(names(count2),collapse=",")))
	} else{
		id_fam<-rbind(id_fam, c(cluster[i,1],NA,NA))
	}	
}
id_fam<-id_fam[,c(2,1,3)]

# edges with a combined weight below the given threshold are drawn transparently 
combined.weight<-network[,"Sum.of.branches"]*network[,"connection.weight"]

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
	#left<-length(grep(paste0(as.character(anno[i,1]),","),cluster[,2]))
	#right<-length(grep(paste0(as.character(anno[i,1]),","),cluster[,3]))
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


# predominant distance from sRNA
left_list<-lapply(cluster[,2], function(x){
		tmp<-strsplit(x,",")[[1]]
		tmp
	}
	)

right_list<-lapply(cluster[,3], function(x){
		tmp<-strsplit(x,",")[[1]]
		tmp
	}
	)	
	
position_integer<-c()
for(i in 1:nrow(anno)){
	temp_left<-unlist(lapply(left_list,function(x){
		return(which(x==anno[i,1]))
	}))
	
	temp_right<-unlist(lapply(right_list,function(x){
		return(which(x==anno[i,1]))
	}))
	
	if(length(temp_left)>=length(temp_right)){
		pos<-table(temp_left)
		pos<-sort(pos, decreasing=T)
		pos<-as.numeric(names(pos)[1])-1
	}	else {
		pos<-table(temp_right)
		pos<-sort(pos, decreasing=T)
		pos<-as.numeric(names(pos)[1])-1
	}	
	if(length(pos)==0){
		pos<-NA
	}
	position_integer<-c(position_integer,pos)
}

anno<-cbind(anno,position_integer)


not_empty<-which(anno[,"filtered_names"]!="")

anno[not_empty,"filtered_names"]<-paste0(anno[not_empty,"filtered_names"],"_", anno[not_empty,"position_integer"])
anno[,"cluster_id"]<-paste0("'",anno[,"cluster_id"],"_",anno[,"position_integer"],"'")

uid<-names(clus_list)

for(i in 1:length(uid)){
	tmp<-rep(NA, ncol(anno))
	tmp[1]<-paste0("sRNA_",i)
	tmp[5]<-max(as.numeric(anno[,"PageRank"]),na.rm=T)*0.6
	tmp[6]<-"sRNA"
	tmp[7]<-"no" 
	anno<-rbind(anno,tmp)
}

write.table(network, file=paste("post_processed_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno, file=paste("post_processed_",anno_file, sep=""),sep="\t", row.names=F, quote=F)


# drawing the detected families as un-connected sub-networks
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
anno2<-c()
for(i in 1:min(length(uid),max_synt=max_synt)){
	tmp<-which(network[,"cluster_id"]==uid[i])
	tmp<-network[tmp,]
	
	if(only_sig_nodes==TRUE){
		tmp<-check_sig(tmp)
	}
	
	notcon<-not_connected2sRNA(tmp)
	if(length(notcon)>0){
		notcon<-notcon[which(tmp[match(notcon,tmp[,1]),"Sum.of.branches"]>=thres_anno)]
		if(length(notcon)>0){
			notcon<-notcon[order(tmp[match(notcon,tmp[,1]),"Sum.of.branches"], decreasing=T)]
			path<-surf(network, notcon)
			if(length(path)>1){
				#find path with highest overlap to sub network
				ov<-unlist(lapply(path, function(x){
					tmp<-length(match(x,c(tmp[,1],tmp[,2])))
				}))
				ov<-which(ov==max(ov))
				path<-path[ov]
				if(length(path)>1){
					# length of path to srna
					le<-unlist(lapply(path, length))
					le<-which(le==min(le))
					path<-path[le]
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
		anno2<-rbind(anno2,tmp_anno)
	}
}

for(i in 1:length(uid)){
	tmp<-rep(NA, ncol(anno2))
	tmp[1]<-paste0("sRNA_",i)
	tmp[5]<-max(as.numeric(anno2[,"PageRank"]),na.rm=T)*0.6
	tmp[6]<-"sRNA"
	tmp[7]<-"no" 
	anno2<-rbind(anno2,tmp)
}


write.table(network2, file=paste("post_processed2_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno2, file=paste("post_processed2_",anno_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(id_fam, file="ids2synt_families.txt",sep="\t", row.names=F, quote=F)





tnet<-which(network2[,"transparent"]=="yes")
if(length(tnet)>0){
	network2<-network2[-tnet,]
}

write.table(network2, file=paste("post_processed_wo_trans_",network_file, sep=""),sep="\t", row.names=F, quote=F)
write.table(anno2, file=paste("post_processed_wo_trans_",anno_file, sep=""),sep="\t", row.names=F, quote=F)









