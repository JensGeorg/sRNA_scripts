#!/usr/bin/Rscript

# dependencies: 
require(stringi)
require(phangorn)
require(ggtree)
#genbank_groper_sqliteDB_ver00.py
#CALL:
#R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r --args write_files=FALSE threads=10  synteny_window=3000 script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney/syntney.db < ~/media/jens@margarita/Syntney/Rfam_db.fasta
#R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r --args write_files=FALSE threads=10 filename=~/media/jens@margarita/Syntney/testfiles/inputForJens.fasta  synteny_window=3000 script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney/new.db


#R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r --args write_files=FALSE threads=10 filename=~/media/jens@margarita/Copra2_paper/Glassgo/RyhB.fa  synteny_window=3000 script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney_db/synt.db
#R --slave -f  ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r --args write_files=FALSE threads=10 filename=~/media/jens@margarita/Syntney/testfiles/Spot42.fa  synteny_window=5000 script_path=~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/Syntney_db/synt.db

#	python3 Syntney.py -i ~/media/jens@margarita/Copra2_paper/Glassgo/Spot42.fa  -o ~/media/jens@margarita/Copra2_paper/Glassgo/Spot42/ -n cys -r off -d ~/Syntney_db/synt.db -c ~/media/jens@margarita/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r  -s ~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py


#filename<-file('stdin', 'r') # result fasta file from GLASSgo
filename<-"~/Copra2_paper/Glassgo/RyhB/RyhB_ref2.fa"
#filename<-"~/media/jens@margarita/Copra2_paper/Glassgo/RyhB/RyhB_ref2.fa"
#filename<-"~/Copra2_paper/Glassgo/Spot42/Spot42_rev.fa"
script_path<-"~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
script_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
#db_path<-"/media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db"
db_path<-"~/synt.db"
#db_path<-"~/Syntney_db/synt.db"
threads<-30
name<-"sRNA"
write_files<-F
rRNA_existence_threshold<-0.6
up<-150
down<-50
wildcard_ids<-c("U00096.3","CP003069.1")
n_orgs_synt<-200
synteny_window<-5000 # number of bases upstream and downstream of the sRNA that were searched for protein coding genes for the synteny analysis
thres_val<-0.0001
args <- commandArgs(trailingOnly = TRUE) 

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }
 
synteny_window<-as.numeric(synteny_window)
threads<-as.numeric(threads)
write_files<-as.logical(write_files)
rRNA_existence_threshold<-as.numeric(rRNA_existence_threshold)
up<-as.numeric(up)
down<-as.numeric(down)
unspli<-grep(",", wildcard_ids)
if(length(unspli)>0){
	wildcard_ids<-strsplit(wildcard_ids,",")[[1]]
}

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

# Execute functions
fasta<-read.delim(filename, header=F, sep="\t")
closeAllConnections()
fasta<-as.character(fasta[,1])
coor<-export_ncRNA_coordinates(fasta)
coor<-remove_overlapping_homologs(coor)


n<-sort(table(paste(coor[,1],coor[,"name"],sep="_")))
n<-cbind(names(n),n)

write.table(n, file="sRNA_per_org.txt", sep="\t")

#16S RNA
orgs<-unique(coor[,1])
command<-paste("python3 ", script_path, " -s ", db_path, " -rRNA ", paste(orgs, collapse=" "))
print(command)
rRNA<-system(command, intern=T)

startp<-grep(">",rRNA)[1]
rRNA<-rRNA[startp:length(rRNA)]

empty<-grep("No 16srRNA found !!!!!!!",rRNA)
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

bed<-c()

for(i in 1:nrow(coor)){
	if(coor[i,2]=="+"){
		tmp<-c(coor[i,1],as.numeric(coor[i,3])-up,as.numeric(coor[i,4])+down, coor[i,2])
	} else {
		tmp<-c(coor[i,1],as.numeric(coor[i,3])-down,as.numeric(coor[i,4])+up, coor[i,2])
	}
	bed<-rbind(bed,tmp)
}
tmpf<-tempfile()
write.table(bed, file=tmpf, quote=F, sep="\t", row.names=F, col.names=F)
command<-paste("python3 ", script_path, " -s ", db_path, " -pdna ", tmpf, sep="")
bed2<-system(command, intern=T)
sta<-grep(">", bed2)[1]
bed2<-bed2[sta:length(bed2)]


combined_seqs<-c()

for(i in 1:nrow(coor)){
	tmp<-c()
	r<-na.omit(which(paste(">",coor[i,1],sep="")==rRNA)[1])
	if(length(r)>0){
		if( rRNA[r+1]!="No 16sRNA sequence found"){
			tmp<-c(tmp, coor[i,"ID"],coor[i,"Full_header"],coor[i,"sequence"],rRNA[r+1])
			if(coor[i,2]=="+"){
				na<-paste(">",coor[i,1],"_",as.numeric(coor[i,3])-up,"_",as.numeric(coor[i,4])+down,sep="")
			}
			if(coor[i,2]=="-"){
				na<-paste(">",coor[i,1],"_",as.numeric(coor[i,3])-down,"_",as.numeric(coor[i,4])+up,sep="")
			}
			pr<-na.omit(which(bed2==na)[1])
			if(length(pr)>0){
				#print(i)
				tmp<-c(tmp,bed2[pr+1])		
				combined_seqs<-rbind(combined_seqs,tmp)
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
align16<-system(paste("mafft --thread 40 --retree 2 --maxiterate 0 --quiet --inputorder ", temp_fas, " > ", temp_fas2,seq=""))

ribo<-read.phyDat(temp_fas2, format="fasta", type="DNA")
dm <- dist.ml(ribo, model="F81")
fitJC<- upgma(dm)
weight16<-tree_weights(fitJC, method="clustal")
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
align_sRNA<-system(paste("mafft --thread 40 --retree 2 --maxiterate 0 --quiet --inputorder ", temp_fas, " > ", temp_fas2,seq=""))
srna<-read.phyDat(temp_fas2, format="fasta", type="DNA")
dm <- dist.ml(srna, model="F81")
fitJC<- upgma(dm)
weight_sRNA<-tree_weights(fitJC, method="clustal")
dist_ham_srna<-dist.hamming(srna)
dist_ham_srna<-as.matrix(dist_ham_srna)

# combined weights
weight_comb<-apply(cbind(weight16,weight_sRNA), 1, function(x){return(exp(mean(log(x))))})
weight_comb<-weight_comb/max(weight_comb)


# alignment of promoters
prom<-c()
for(i in 1:nrow(combined_seqs)){
	prom<-c(prom,paste(">", combined_seqs[i,1], sep=""))
	prom<-c(prom, combined_seqs[i,5])
}

temp_fas<-tempfile()
temp_fas2<-tempfile()
writeLines(prom,con=temp_fas)
align_prom<-system(paste("mafft --thread 40 --retree 2 --maxiterate 0 --quiet --inputorder ", temp_fas, seq=""),intern=TRUE)

header<-grep(">",align_prom)
seq4<-vector("list", length(header))
names(seq4)<-gsub(">","",align_prom[header])
for(j in 1:(length(header)-1)){
	tmp<-align_prom[(header[j]+1):(header[j+1]-1)]
	tmp<-paste(tmp,collapse="")
	seq4[[j]]<-strsplit(tmp,"")[[1]]
}
tmp<-align_prom[(header[length(header)]+1):length(align_prom)]
tmp<-paste(tmp,collapse="")
seq4[[length(seq4)]]<-strsplit(tmp,"")[[1]]


seq4<-do.call(rbind,seq4)
seq4<-toupper(seq4)
weight_comb<-weight_comb*nrow(combined_seqs)
new_num<-sum(weight_comb)
score<-c()
m_score<-c()
mat<-matrix(,ncol(seq4),4)
colnames(mat)<-c("A","T","C","G")
rownames(mat)<-1:ncol(seq4)
for(j in 1:ncol(seq4)){
	mat[j,1]<-sum(weight_comb[which(seq4[,j]=="A")])#/new_num2
	mat[j,2]<-sum(weight_comb[which(seq4[,j]=="T")])##/new_num2
	mat[j,3]<-sum(weight_comb[which(seq4[,j]=="C")])#/new_num2
	mat[j,4]<-sum(weight_comb[which(seq4[,j]=="G")])#/new_num2
}
		
mat<-capture.output(write.table(mat, file=stdout(), sep=" ", quote=F,row.names = T ))
mat[1]<-paste("P0 ", mat[1],sep="")
mat<-c("ID any_old_name_for_motif_1","BF species_name_for_motif_1",mat,"XX")
writeLines(mat, con="MAt.pssm")
writeLines(align_prom, con="align.fa")
id_fam<-read.csv("ids2synt_families.txt", sep="\t", header=T)

po<-match(id_fam[,2], coor[,"ID"])
na<-which(is.na(po)==F)
po<-na.omit(po)
fam<-matrix(,nrow(coor),3)
fam[po,1]<-id_fam[na,1]
fam[po,2]<-id_fam[na,2]
fam[po,3]<-id_fam[na,3]


co3<-cbind(coor,fam)
write.table(co3,file="annotated_homologs.txt", sep="\t")

na<-gsub(" .*","",coor[,"name"])
nau<-unique(na)
numb_list<-vector("list",length(nau))
names(numb_list)<-nau
for(i in nau){
	tmp<-which(na==i)
	ac<-unique(coor[tmp,1])
	numb<-c()
	for(j in ac){
		numb<-c(numb,length(which(coor[,1]==j)))
	}
	numb_list[[i]]<-numb
}
################


# exclude_similars<-function(dis, thres=0.01,ooi){
# i<-1


# for(j in 1:length(ooi)){
	# o<-match(ooi, gsub("\\..*","",rownames(dis)))
	# nam_ref<-colnames(dis)[o]
	# temp<-which(dis[,o[j]]<=thres)
	
	# temp<-setdiff(temp,o)
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		
	# }
	# }



# while(nrow(dis)>i){
	# o<-match(ooi, gsub("\\..*","",rownames(dis)))
	# nam<-colnames(dis)[i]
	# o2<-match(nam,colnames(dis))
	# temp<-which(dis[,i]<=thres)
	# temp<-setdiff(temp,c(o2,o))
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		
	# }
	# i<-i+1
# }
# dis
# }


exclude_similars<-function(dis,dis2, thres=0.01, thres2=0.01,ooi, id_fam=NULL){
i<-1

	if(is.null(id_fam)==F){
		ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
		na<-which(is.na(ids))
		if(length(na)>0){
			dis<-dis[-na,-na]
			dis2<-dis2[-na,-na]
		}
		
		ooi2<-c()
		for(jj in ooi){
			ooi2<-c(ooi2,grep(jj, id_fam[,2]))
		}
		ooi<- id_fam[ooi2,2]
	}
	

	for(j in 1:length(ooi)){
		o<-c()
		for(jj in ooi){
			o<-c(o,which(gsub("\\..*","",rownames(dis))==gsub("\\..*","",jj)))
		}
		#match(ooi, gsub("\\..*","",rownames(dis)))
		#nam_ref<-colnames(dis)[o]
		
		if(is.null(id_fam)){
			temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 )
			temp<-setdiff(temp,o)
		} else {
			id_ref<-id_fam[which(id_fam[,2]==ooi[i]),1]
			#if(length(id_ref)>1){print(c(id_ref,o))}
			ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
			temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 & ids==id_ref[1] )
			temp<-setdiff(temp,o)
		}
		if(length(temp)>0){
			dis<-dis[-temp,-temp]
			dis2<-dis2[-temp,-temp]
		}
		
	}


while(nrow(dis)>i){
	oo<-c()
	for(jj in ooi){
		oo<-c(oo,which(gsub("\\..*","",rownames(dis))==gsub("\\..*","",jj)))
	}
	nam<-colnames(dis)[i]
	o<-grep(nam, rownames(dis))
	if(is.null(id_fam)){
		temp<-which(dis[,i]<=thres & dis2[,i]<=thres2)
		temp<-setdiff(temp,c(oo,i))
	}	else {
		id_ref<-id_fam[which(id_fam[,2]==nam),1]
		#if(length(id_ref)>1){print(c(id_ref,nam))}
		ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
		temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 & ids==id_ref[1] )
		temp<-setdiff(temp,c(oo,i))
	}
	
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		dis2<-dis2[-temp,-temp]
	}
	i<-i+1
}
dis
}






# ########################
# # function to exclude very similar organism based on a phylogentic tree to reduce complexity
# exclude_similars<-function(dis,dis2, thres=0.01, thres2=0.01,ooi, id_fam=NULL){
# i<-1

# if(is.null(id_fam)==F){
	# ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
	# na<-which(is.na(ids))
	# if(length(na)>0){
		# dis<-dis[-na,-na]
		# dis2<-dis2[-na,-na]
	# }
# }

# o<-grep(ooi, rownames(dis))
# nam<-colnames(dis)[o]
	# if(is.null(id_fam)){
		# temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 )
	# } else {
		# id_ref<-id_fam[which(id_fam[,2]==ooi),1]
		# #if(length(id_ref)>1){print(c(id_ref,o))}
		# ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
		# temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 & ids==id_ref[1] )
	# }
	# nam<-na.omit(match(nam, rownames(dis)))
	# nam<-na.omit(match(nam, temp))
	# if(length(nam)>0){
		# temp<-temp[-nam]
	# }
	
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		# dis2<-dis2[-temp,-temp]
	# }
	



# while(nrow(dis)>i){
	
	# nam<-colnames(dis)[i]
	# o<-grep(nam, rownames(dis))
	# if(is.null(id_fam)){
		# temp<-which(dis[,i]<=thres & dis2[,i]<=thres2)
	# }	else {
		# id_ref<-id_fam[which(id_fam[,2]==nam),1]
		# #if(length(id_ref)>1){print(c(id_ref,nam))}
		# ids<-id_fam[match(rownames(dis),id_fam[,2]),1]
		# temp<-which(dis[,o]<=thres & dis2[,o]<=thres2 & ids==id_ref[1] )
	# }
	# nam<-na.omit(match(nam, rownames(dis)))
	# nam<-na.omit(match(nam, temp))
	# if(length(nam)>0){
		# temp<-temp[-nam]
	# }
	
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		# dis2<-dis2[-temp,-temp]
	# }
	# i<-i+1
# }
# dis
# }

steps<-0
thres_vect<-c()
a<-dist_ham_srna
dim_vect<-c()
if(dim(a)[1] > n_orgs_synt){
	while(is.null(a) | dim(a)[1] > n_orgs_synt){
		a<-exclude_similars(dist_ham_srna,dist_ham_ribo ,thres=thres_val, thres2=thres_val, id_fam=id_fam,ooi=wildcard_ids)
		thres_vect<-c(thres_vect,thres_val)
		dims<-dim(a)[1]
		names(dims)<-thres_val
		dim_vect<-c(dim_vect,dims)
		thres_val<-thres_val*10		
	}
	while(steps < 10 & dim(a)[1] != n_orgs_synt){
		thres_vect<-thres_vect[(length(thres_vect)-1):length(thres_vect)]
		if(dim(a)[1]<n_orgs_synt){
			thres_val<-thres_vect[2]-abs(thres_vect[1]-thres_vect[2])/2
		} else {
			thres_val<-thres_vect[2]+abs(thres_vect[1]-thres_vect[2])/2
		}
		thres_vect<-c(thres_vect,thres_val)
		a<-exclude_similars(dist_ham_srna,dist_ham_ribo ,thres=thres_val, thres2=thres_val, id_fam=id_fam,ooi=wildcard_ids)
		steps<-steps+1
		dims<-dim(a)[1]
		names(dims)<-thres_val
		dim_vect<-c(dim_vect,dims)
	}
	thres_val<-as.numeric(names(which(abs(dim_vect-n_orgs_synt)==min(abs(dim_vect-n_orgs_synt)))[1]))
	a<-exclude_similars(dist_ham_srna,dist_ham_ribo ,thres=thres_val, thres2=thres_val, id_fam=id_fam,ooi=wildcard_ids)
}






c2<-combined_seqs[match(colnames(a),combined_seqs[,"ID"]),]

prom2<-c()
for(i in 1:nrow(c2)){
	prom2<-c(prom2,paste(">", c2[i,1], sep=""))
	prom2<-c(prom2, c2[i,5])
}

temp_fas<-tempfile()
temp_fas2<-tempfile()
writeLines(prom2,con=temp_fas)
align_prom2<-system(paste("mafft --thread 40 --maxiterate 1000 --localpair --quiet --inputorder  ", temp_fas," > ", temp_fas2,seq=""),intern=TRUE)




pr<-read.phyDat(temp_fas2, format="fasta", type="DNA")
dm <- dist.ml(pr, model="F81")




treeNJ <- NJ(dm)
mt <- modelTest(pr, tree=treeNJ, multicore=TRUE)
mt[order(mt$AICc),]
bestmodel <- mt$Model[which.min(mt$AIC)]
env = attr(mt, "env")
fitStart = eval(get(bestmodel, env), env)
fitJC = optim.pml(fitStart, rearrangement = "stochastic",optGamma=TRUE, optInv=TRUE, model="GTR")


 # fitStart = pml(treeNJ, pr, k=4)
 # fitJC = optim.pml(fitStart, model="GTR", optGamma=T, rearrangement="stochastic",ratchet.par = list(iter = 5L, maxit = 20L, prop = 1/3),control = pml.control(epsilon = 1e-08, maxit = 10,
 # trace = 1L))

tree<-fitJC$tree
#tree<-fit$tree
nam<-match(tree$tip.label, coor[,"ID"])
nam_id<-match(tree$tip.label, id_fam[,2])
tree$tip.label<-paste(coor[nam,"name"],coor[nam,"ID"],"|",id_fam[nam_id,1],";",id_fam[nam_id,3],sep="")
lab2<-id_fam[nam_id,1]


# require(ggtree)
# pdf("tree_test.pdf")
# ggtree(tree)  + geom_tiplab(size=0.8)
# dev.off()


############ tree collapsing

#p<-ggtree(tree)  + geom_tiplab(size=0.8) #+ geom_cladelabel(node=333, label='Shewanella', angle=-95, hjust=.5, fontsize=8)

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
	tips<-out[which(is.element(out,topo[,1])==F)]
	nodes<-setdiff(out,tips)
	return(list(tips,nodes))	
}


root<-setdiff(tree[[1]][,1],tree[[1]][,2])
set<-root
tips<-setdiff(tree[[1]][,2],tree[[1]][,1])
tips2<-tips
visited_nodes<-root
out_nodes<-c()

while(length(set)>0){
	set_temp<-c()
	for(i in 1:length(set)){
		la<-leafs(set[i],tree)[[1]]
		no<-leafs(set[i],tree)[[2]]
		#print(c(set[i],unique(lab2[la])))
		if(length(unique(lab2[la]))==1){
			out_nodes<-c(out_nodes,set[i])
			#print(c(set[i],unique(lab2[la])))
			#tips<-setdiff(tips,la)
			#la<-setdiff(la,tips)
			visited_nodes<-c(visited_nodes,no,la)
		} else {
			pos1<-na.omit(tree[[1]][which(tree[[1]][,1]==set[i]),2])
			set_temp<-c(set_temp,pos1)	
			set_temp<-setdiff(set_temp,visited_nodes)
		}
	}
	set<-set_temp
}
out_nodes<-setdiff(out_nodes,tips2)
color<-c("#117893","#D13B40","#7AB800","#477831","#EAAC9D","#FFAF12","#4EC5A5","#666666","#A2553A","#72617D")

if(length(color)<length(out_nodes)){
	nu<-length(out_nodes)-length(color)
	color2<-c(color, rep(1,nu))
} else{
	color2<-color
}

st2<-"p<-ggtree(tree, layout='circular') "#+ geom_tiplab(size=0.8, col=color2[cl] )"# + geom_tiplab(size=0.8, col=color2[cl])"#, layout='circular') "#  "
st<-"p"
cl<-c()
ti<-c()
na<-c()
an<-c()
for(i in 1:length(out_nodes)){
	la<-leafs(out_nodes[i],tree)[[1]]
	cl<-c(cl,lab2[la])
	ti<-c(ti,la)
	#print(c(out_nodes[i], unique(lab2[la])))
	st<-paste(st, " %>% collapse(node=c(", out_nodes[i], "),color='black',fill='",color2[as.numeric(unique(lab2[la]))],"', 'max') ", sep="")
	tmp<-unique(gsub(" .*","",tree$tip.label[la]))
	if(length(tmp)>0){
		tmp<-paste(sort(tmp), collapse="\n")
	}
	an<-c(an,tmp)
	if(is.element(out_nodes[i],tips2)==F){
		st2<-paste(st2, " + geom_cladelabel(node=",out_nodes[i],",label='",(tmp),"', align=T, fontsize=1,color='",color2[as.numeric(unique(lab2[la]))],"')", sep="")
	}
}


cl<-cl[order(ti)]
 eval(parse(text=st2))

#p<-ggtree(tree)  + geom_tiplab(size=0.8, col=color2[cl]) #+ geom_cladelabel(node=333, label='Shewanella', angle=-95, hjust=.5, fontsize=8)

 #eval(parse(text=st))

pdf("sRNA_circ.pdf")

 eval(parse(text=st))

dev.off()


st2<-"p<-ggtree(tree) + geom_tiplab(size=0.8)"#, col=color2[cl] )"# + geom_tiplab(size=0.8, col=color2[cl])"#, layout='circular') "#  "
st<-"p"
cl<-c()
ti<-c()
na<-c()
an<-c()
for(i in 1:length(out_nodes)){
	la<-leafs(out_nodes[i],tree)[[1]]
	cl<-c(cl,lab2[la])
	ti<-c(ti,la)
	#print(c(out_nodes[i], unique(lab2[la])))
	st<-paste(st, " %>% collapse(node=c(", out_nodes[i], "),fill='",color2[as.numeric(unique(lab2[la]))],"', 'max') ", sep="")
	tmp<-unique(gsub(" .*","",tree$tip.label[la]))
	if(length(tmp)>0){
		tmp<-paste(sort(tmp), collapse="\n")
	}
	an<-c(an,tmp)
	if(is.element(out_nodes[i],tips2)==F){
		st2<-paste(st2, " + geom_cladelabel(node=",out_nodes[i],",label='",(tmp),"', align=T, fontsize=1,color='",color2[as.numeric(unique(lab2[la]))],"')", sep="")
	}
}


cl<-cl[order(ti)]
 eval(parse(text=st2))


pdf("sRNA_reference.pdf")

p

dev.off()





# pdf("tree_circ.pdf")
# p<-ggtree(tree)  + geom_tiplab(size=0.8) #+ geom_cladelabel(node=333, label='Shewanella', angle=-95, hjust=.5, fontsize=8)

 # p %>% collapse(node=c(317), "mixed") %>% collapse(313, 'mixed') %>% collapse(292, 'mixed')%>%
# collapse(301, 'mixed') %>% collapse(304, 'mixed') %>% collapse(257, 'mixed') %>%
# collapse(262, 'mixed') %>% collapse(281, 'mixed') %>% collapse(333, 'mixed') %>%
# collapse(345, 'mixed') %>% collapse(365, 'mixed') %>% collapse(360, 'mixed')  %>%
# collapse(354, 'mixed')  %>% collapse(371, 'mixed')  %>%
# collapse(418, 'mixed') %>% collapse(395, 'mixed') %>% collapse(391, 'mixed')  %>%
# collapse(411, 'mixed') %>%  collapse(405, 'mixed') %>%
# collapse(439, 'mixed') %>% collapse(229, 'mixed') %>% collapse(359, 'mixed') %>%
# collapse(348, 'mixed') %>% collapse(404, 'mixed') %>% collapse(402, 'mixed') %>%
# collapse(431, 'mixed') %>% collapse(419, 'mixed') %>% collapse(366, 'mixed') %>%
# collapse(416, 'mixed') %>% collapse(436, 'mixed') %>% collapse(241, 'mixed') %>%
# collapse(430, 'mixed') %>% collapse(226, 'mixed') %>% collapse(425, 'mixed') %>%
# collapse(327, 'mixed') %>% collapse(260, 'mixed') %>% collapse(276, 'mixed') %>%
# collapse(249, 'mixed') %>% collapse(426, 'mixed')  %>%    collapse(352, 'mixed') %>% collapse(353, 'mixed') 
# dev.off()


# pdf("tree_with_nodes_cl2.pdf")
# ggtree(tree)  + geom_tiplab(size=0.8) + geom_label2(aes(subset=!isTip, label=node), col=2,size=1,label.padding=unit(0, "lines"),label.size = 0)
# dev.off()

# save.image()

coor2<-match(fitJC$tree$tip.label, coor[,"ID"])
coor2<-coor[coor2,]
colnames(id_fam)<-c("synteny_cluster","ID","all_cluster")
pos<-match(coor2[,"ID"],id_fam[,2])
coor2<-cbind(coor2,id_fam[pos,])
coor2<-coor2[order(coor2[,"synteny_cluster"]),]
s<-as.numeric(coor2[,3])
e<-as.numeric(coor2[,4])
m<-round(s+(e-s)/2,digits=0)
synteny_window<-5000
wi<-rep(synteny_window,nrow(coor2))

#coor3<-cbind(paste(coor[,1],"_",coor[,3],sep=""),coor[,1],m,wi)
coor3<-cbind(coor2[,"ID"],coor2[,1],m,wi)

coordinates<-tempfile()
write.table(coor3,file=coordinates, sep="\t", row.names=F, col.names=F, quote=F)

command<-paste("python3 ", script_path, " -s ", db_path, " -a ", coordinates)
print(command)
dat<-system(command, intern=T)
empty<-which(dat=="")
if(length(empty)>0){
	dat<-dat[-empty]
}
dat<-do.call(rbind,strsplit(dat,"\t")) 

na<-grep("no annotation", dat[,3])
na2<-na
if(length(na)>0){
	na<-dat[na,1]
}

na<-c(setdiff(coor3[,1],unique(dat[,1])),na)
if(length(na)>0){
	tmp<-unlist(lapply(na, function(x){ return(which(coor[,"ID"]==x))}))
	coor2<-coor2[-tmp,]
	dat<-dat[-na2,]
}
unlink(coordinates)

ids<-unique(dat[,1])
id2<-paste(coor2[,"ID"])


ids<-na.omit(ids[match(id2,ids)])

out<-vector("list", length(ids))
names(out)<-ids

for(i in 1:length(ids)){
	tmp<-which(dat[,1]==ids[i])
	srna<-match(ids[i],id2)
	if(is.na(srna)==FALSE){
		stra<--1
		if(coor2[srna,2]=="+"){
		  stra<-1
		}
		temp_out<-cbind(dat[tmp,7],dat[tmp,5],dat[tmp,6],dat[tmp,3],dat[tmp,4],dat[tmp,8])
		colnames(temp_out)<-c("strand","start","end","gene_name","locus_tag","AA_sequence")
		ma<-max(as.numeric(temp_out[,3]))
		mi<-min(as.numeric(temp_out[,2]))
		aa<-as.numeric(temp_out[,2])-mi
		bb<-as.numeric(temp_out[,3])-mi
		s_srna<-min(as.numeric(coor2[srna,3:4]))-mi
		e_srna<-max(as.numeric(coor2[srna,3:4]))-mi
		nan<-which(is.na(temp_out[,"locus_tag"]))
		na<-which(temp_out[,"locus_tag"]=="na")
		nan<-unique(c(na,nan))
		if(length(nan)>0){
			temp_out[nan,"locus_tag"]<-unlist(lapply(length(nan),rand_extension,Accession=ids[i]))
		}
		temp_out<-data.frame(temp_out,aa,bb,rep(s_srna,nrow(temp_out)),rep(e_srna,nrow(temp_out)),rep(stra,nrow(temp_out)),rep(coor2[srna,5],nrow(temp_out)))
		out[[i]]<-temp_out
	}
}
filen<-tempfile()
tagtable<-locus_tag2org(out)
get_prot_fasta3(out, filen)
cd<-cdhit_run(fasta=filen, psi=F,thres=0.4, threads=threads)
cd<-proc_cdhit(cd)


plot_function4<-function(out, cdhit_result, wind=3000, outformat="cairo_pdf", fasta="sRNA"){ 
  cdl<-unlist(lapply(cdhit_result,length))
  one<-which(cdl==1)
  more<-which(cdl>1)
  nl<-length(more)
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  while(length(more)>length(color)){
	color<-c(color,color)
  }
  if(length(more)>0){
    collist<- sample(color, nl)
    for(i in 1:length(more)){
      names(cdhit_result)[more[i]]<-collist[i]
    }
  }
  names(cdhit_result)[one]<-"grey"
  numb<-ceiling(length(out)/5)
  
  
  
  le<-length(out)
  
   count<-4*le+1
	if(outformat=="pdf"){
	nam<-paste(fasta,"_synteny.pdf",sep="")
	pdf(nam,width=12, height=le*0.9,  useDingbats=F)
	#cairo_pdf(nam,width=11,  height = 1.1*le)
  }
  if(outformat=="png"){
	nam<-paste(fasta,"_synteny.png",sep="")
	png(nam,width=11,  height = 1.1*le,units="in", res=100)
  }
  if(outformat=="svg"){
	nam<-paste(fasta,"_synteny.svg",sep="")
	svg(nam,width=11,  height = 1.1*le)
  }
	for( i in 1:le){
		if(out[[i]][1,11]==-1){
		temp<-out[[i]]
		
		temp[,11]<-rep(1,nrow(temp))
		
		ma<-max(temp[,8])
		plus<-which(as.character(temp[,1])=="+")
		minus<-which(as.character(temp[,1])=="-")
		ord<-rep(NA,nrow(temp))
		levels(temp[,1])<-c("+","-")
		if(length(plus)>0){
			temp[plus,1]<-"-"
			ord[plus]<-"-"
		}
		if(length(minus)>0){
			temp[minus,1]<-"+"
			ord[minus]<-"+"
		}
			for(j in 1:nrow(temp)){
				temp[j,8]<-ma-out[[i]][j,7]
				temp[j,7]<-ma-out[[i]][j,8]
				temp[j,10]<-ma-out[[i]][j,9]
				temp[j,9]<-ma-out[[i]][j,10]
				}
			out[[i]]<-temp	
		}
		
	}
	
	#mas<-out[[1]][1,9]
	mas<-wind
	 for(i in 1:le){
		 te<-out[[i]]
		 s2<-mas-as.numeric(te[1,9])
		for(ii in 1:nrow(te)){
        
         
		 #print(as.numeric(te[1,9])[1])
        #s<-c(s,s2+mas)
         #se<-c(se,s2+mi)
         te[ii,7]<-as.numeric(te[ii,7])+s2
         te[ii,8]<-as.numeric(te[ii,8])+s2
         te[ii,9]<-as.numeric(te[ii,9])+s2
         te[ii,10]<-as.numeric(te[ii,10])+s2
         
		}
		out[[i]]<-te
     }
	
  #
  
  
  
  d_vect<-c()
  s_vect<-c()
  
	plot(1,1, type="n",xlim=c(1,wind*2),ylim=c(1,count),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
    for(i in 1:le){
    
        temp<-out[[i]]
		nam3<-names(out)[i]
        mi<-min(as.numeric(temp[,7]))
        ma<-max(as.numeric(temp[,8]))

        lines(c(mi,ma),c(1+count-i*4,1+count-i*4))
        lines(c(mi,ma),c(2+count-i*4,2+count-i*4))
        for(j in 1:nrow(temp)){

          m<-(as.numeric(temp[j,8])-as.numeric(temp[j,7]))/2+as.numeric(temp[j,7])
          mm<-as.numeric(temp[j,8])-as.numeric(temp[j,7])
          n<--1
          if(temp[j,1]=="+"){
            n<-1
          }
          if(n==-1){
            n<-0
          }
          color<-"white"
          tcolor<-na.omit(grep(temp[j,5],cdhit_result))
          tcolor2<-na.omit(grep(temp[j,4],cdhit_result))
          tcolor<-c(tcolor,tcolor2)
          if(length(tcolor)>0){
             color<-names(cdhit_result)[tcolor]
           }
          rect(as.numeric(temp[j,7]),0.5+n+count-i*4,as.numeric(temp[j,8]),1.5+n+count-i*4, col=color)
          text(m,n+count-i*4+1+0.25,gsub("\"","",gsub("NA ","",gsub("orf.*","",temp[j,4]))),cex=0.4,font=4)
		  text(m,n+count-i*4+1-0.25,gsub("\"","",temp[j,5]),cex=0.4,font=4)
        }
        n<-as.numeric(temp[1,11])
        if(n==-1){
          n<-0
        }
        rect(as.numeric(temp[j,9]),0.5+n+count-i*4,as.numeric(temp[j,10]),1.5+n+count-i*4, col=2)
        text(wind,count-i*4,paste(nam3, temp[1,12], sep="_"), cex=0.8, font=4)
  	
    }
 dev.off()	

}

plot_function4(out,cd,wind=synteny_window,outformat="pdf", fasta="Spot42")


########## Copra selection



