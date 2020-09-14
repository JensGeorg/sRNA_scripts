
# GLASSgo postprocessing script

# dependencies: 
## CDhit
require(seqinr)
#require(RSQLite)
#require(rentrez)
require(phangorn)
#require(stringi)
# parameters from function call:


#R --slave -f /media/cyano_share/data/WM/sRNA/coprarna_selection.r --args filename=esrA   name=EsrA_15 ooi=NC_002516.2 maxorgs=15 maxdis=1.5 clustervalue=1 mindis=0.001

#CALL:
#R --slave -f  /home/jens/jensSicherung/GLASSgo2/GLASSgo_postprocessing_10_steffen_databank_sqlite.r --args filename=debug.fasta duplicates_allowed=TRUE synteny_window=3000 refpath=/media/cyano_share/data/GLASSgo_postprocessing/accession_to_refseq cop_path=/media/cyano_share/data/GLASSgo_postprocessing/CopraRNA_available_organisms.txt  name=testfile coprarna_compatible=TRUE
#R --slave -f  /home/jens/jensSicherung/GLASSgo2/GLASSgo_postprocessing_6.r --args filename="sRNA.txt" duplicates_allowed=FALSE synteny_window=3000  name=FnrS coprarna_compatible=TRUE ooi=NC_000913



#filename<-file('stdin', 'r') # result fasta file from GLASSgo
#filename<-"~/media/jens@margarita/Syntney/testfiles/error2.fasta"
filename<-"~/Copra2_paper/Glassgo/RyhB/RyhB_ref2.fa"

script_path<-"~/media/jens@margarita/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB_ver01.py"
script_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
#db_path<-"/media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db"
cop_path<-"~/CopraRNA-git/update_kegg2refseq/run/CopraRNA_available_organisms.txt"
db_path<-"~/synt.db"
threads<-30
name<-"sRNA"
write_files<-F
rRNA_existence_threshold<-0.6

synteny_window<-3000 # number of bases upstream and downstream of the sRNA that were searched for protein coding genes for the synteny analysis



#ooi=c("NC_000964.3","NZ_CP018205","NC_007795")
#duplicates_allowed<-F  # if FALSE only one homolog from one organism is plotted
name<-"RsaA"  # name of the investigated sRNA

ooi<-"NC_000913" # organism of interest - important for the CopraRNA organism selection

coprarna_compatible<-T # organsims are filtered based on the CopraRNA_available_organisms.txt file

maxorgs<-15 # number of orgs for CopraRNA prediction
mindis<-0.001 # organisms which have a lower cophonetic distances to each other in a 16S phylogentic tree are excluded to reduce complexity. A respective reference organism is kept
maxdis<-1.5 # # organisms which have a higher cophonetic distances to each the ooi are excluded.
closeorgs<-3 # number of closest relatives to the ooi for CopraRNA
clustervalue<-2 # the higher the higher the resolution of tree sub-groups
# S. oneidensis, V. cholerae, Aeromonas hydrofila, Y. pseudotubecolosis YPIII, S. enterica, Xenorhabdus Nematophila,  Aliivibrio salmonicida
wildcard<-c("NC_004347","NC_002505","NZ_CP006870","NC_010465","NC_016810","NC_014228","NC_011312")
wildcard<-c("NC_004347","NC_002505","NZ_CP006870","NC_010465","NC_016810","NC_014228")
refseq_required=T
outorg="NC_000964.3" # genome accssion of outgroup organsim for phylogenetic tree, if "internal" organism with highest distance to ooi is selected as outgroup
min_size<-2  # minimal size of tree-based cluster
mixed_sources=T # True if several GLAssgo outputs are pasted together
reduced_svg=F
random_extension=F # make locus_tags unique by adding random extensions
full_tree=F
exact_tree=T
rRNA_phylogeny<-T

args <- commandArgs(trailingOnly = TRUE) 

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }
 full_tree<-as.logical(full_tree)
duplicates_allowed<-as.logical(duplicates_allowed)
synteny_window<-as.numeric(synteny_window)
mindis<-as.numeric(mindis)
maxorgs<-as.numeric(maxorgs)
maxdis<-as.numeric(maxdis)
clustervalue<-as.numeric(clustervalue)
closeorgs<-as.numeric(closeorgs)
exact_tree<-as.logical(exact_tree)
rRNA_phylogeny<-as.logical(rRNA_phylogeny)


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
	b<-tmp2[2]
	a<-strsplit(a,"c")[[1]]
	strand<-"+"
	if(length(a)==2){
		st<-b
		en<-a[2]
		strand<-"-"
	}else{
		st<-a
		en<-b
	}
	out<-c(id, strand, st,en,name)
	out
}


export_ncRNA_coordinates<-function(x){ 
	header_row <- grep(">", x)
	headers <- as.character(x[header_row])
	seqs<-c()
	for(i in 1:(length(header_row)-1)){
		s<-header_row[i]+1
		e<-header_row[i+1]-1
		seqs<-c(seqs,paste(x[s:e],collapse=""))
	}
	seqs<-c(seqs,paste(x[(header_row[length(header_row)]+1):length(x)],collapse=""))
	head1<-headers[1]
	seq1<-seqs[1]
	headers<-headers[2:length(headers)]
	#seqs<-as.character(x[header_row+1])
	seqs<-seqs[2:length(seqs)]
	tmp<-do.call(rbind,lapply(headers,split_glassgo))
	tmp<-cbind(tmp,headers,seqs)
	colnames(tmp)<-c("Accesion_number", "Strand","start","end","name","Full_header","sequence")
	out<-list(tmp,c(head1,seq1))
	out
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

to_refseq2<-function(result, refpath="../accession_to_refseq"){
		load(refpath)
		result<-gsub("\\..*","",result)
		temp<-c()
		
		for(i in 1:length(result)){
		temp1<-grep(result[i],ref[,"full_genome_entry"])
		temp<-c(temp,temp1[1])
		}
		temp<-ref[temp,"Chromosomes.RefSeq"]
		refseq<-temp
		refseq
	}

rand_extension<-function(x, Accession){
	ra<-stri_rand_strings(x,length=4,  pattern = "[A-Za-z0-9]")
	temp<-paste(Accession,ra,sep="_")
	#temp<-gsub("\"","",temp)
	
	temp
}


divide<-function(number,candlist,p4,max_step=5,ooi){
	ooi_p<-grep(ooi, colnames(p4))
	#p5<-names(sort(p4[,ooi_p))
	more<-0
	stepsize<-floor(length(candlist)/number)
	stepsize<-min(max_step,stepsize)
	if(number>1){
		sel<-seq(1,length(candlist)*2,by=stepsize)
		sel<-sel[1:number]
		sel<-na.omit(sel)
		more<-number-length(sel)
		if(length(candlist)>1){
			cands<-names(sort(p4[candlist,ooi_p]))
		}
		if(length(candlist)==1){
			cands<-rownames(p4)[candlist]
		}
	}
	if(number==1){
		sel<-1
		#cands<-rownames(p4)[candlist]
		if(length(candlist)>1){
			cands<-names(sort(p4[candlist,ooi_p]))
		}
		if(length(candlist)==1){
			cands<-rownames(p4)[candlist]
		}
	}
	
	sel<-cands[sel]
	if(more>0){
		r<-setdiff(cands,sel)
		sel<-c(sel,r[1:more])
		}	
	#print(c(i,sel))
	sel	
}

# get_prot_fasta3<-function(out){
  # fasta<-c()
  # for(i in 1:length(out)){
    # for(j in 1:nrow(out[[i]])){
        # if(is.na(out[[i]][j,6])==F){
        # temp<-as.character(out[[i]][j,6])
		# na<-as.character(out[[i]][j,5])
		# na<-gsub("\\\"","",na)
         # na<-paste(">",na,sep="")
		 # temp<-c(na,temp)
         # fasta<-c(fasta,temp)
        # }
      # }
    # }
  # write.table(fasta, file="protein_fasta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
# }


# # identify homologous proteins using CDhit
# cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T){
  # wd<-getwd()
  # di<-paste(wd,"/", "psi_out",sep="")
  # dir.create(di)
  # if(psi==T){
    # inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -d 50 -o ", di,"/",outname, " -c ", thres, sep="")
  # }
  # if(psi==F){
    # inp<-paste("cd-hit -i ", fasta,  " -d 50 -o ", di,"/",outname, " -c ",  thres ," -n 2", " -aL 0.6", sep="")
  # }
  # print(inp)
  # system(inp)
  # cd<-paste(di, "/", outname, ".clstr", sep="")
  # cd<-read.delim(cd, header=F, sep="?")
  # cd<-as.character(cd[,1])
  # cd<-gsub("\t"," ", cd)
  # cd
# }

# proc_cdhit<-function(x){ 
  # clustlist<-list()
  # numb<-grep(">Cluster", x)
  # for(i in 1:length(numb)){
    # if(i<length(numb)){
      # end<-numb[i+1]-1
    # }
    # if(i==length(numb)){
      # end<-length(x)
    # }
    # temp<-x[(numb[i]+1):end]
    # temp<-gsub(".*aa, >","",temp)
    # temp<-gsub("\\.\\.\\..*","",temp)
    # clustlist[[i]]<-temp
  # }
  # clustlist	
# }

# plot a pdf visualizing the synteny of predicted sRNA homologs


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

# # identify homologous proteins using CDhit
# cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T, threads=2){
	# tempf<-tempfile()
	# wd<-getwd()
	# di<-paste(wd,"/", "psi_out",sep="")
	# dir.create(di)
	# if(psi==T){
		# inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -d 50 -o ", tempf, " -c ", thres, sep="")
	 # }
	# if(psi==F){
		# inp<-paste("cd-hit -i ", fasta,  " -d 50 -o ",tempf, " -c ",  thres ," -n 2", " -aL 0.6", " -T ", threads, sep="")
	# }
	 # # print(inp)
	 # system(inp, intern=T)
	 # cd<-paste(tempf, ".clstr", sep="")
	 # cd<-readLines(cd)
	 # cd<-as.character(cd)
	 # cd<-gsub("\t"," ", cd)
	 # cd
# }

# proc_cdhit<-function(x){ 
  # clustlist<-list()
  # numb<-grep(">Cluster", x)
  # for(i in 1:length(numb)){
    # if(i<length(numb)){
      # end<-numb[i+1]-1
    # }
    # if(i==length(numb)){
      # end<-length(x)
    # }
    # temp<-x[(numb[i]+1):end]
    # temp<-gsub(".*aa, >","",temp)
    # temp<-gsub("\\.\\.\\..*","",temp)
    # clustlist[[i]]<-temp
  # }
  # clustlist	
# }


# Execute functions



  # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  # while(length(more)>length(color)){
	# color<-c(color,color)
  # }
  # if(length(more)>0){
    # collist<- sample(color, nl)
    # for(i in 1:length(more)){
      # names(cdhit_result)[more[i]]<-collist[i]
    # }
  # }
  # names(cdhit_result)[one]<-"grey"
  # numb<-ceiling(length(out)/5)
  
  
  
  # le<-length(out)
  
   # count<-4*le+1
	# if(outformat=="pdf"){
	# nam<-paste(fasta,"_synteny.pdf",sep="")
	# pdf(nam,width=12, height=le*0.9,  useDingbats=F)
	# #cairo_pdf(nam,width=11,  height = 1.1*le)
  # }
  # if(outformat=="png"){
	# nam<-paste(fasta,"_synteny.png",sep="")
	# png(nam,width=11,  height = 1.1*le,units="in", res=100)
  # }
  # if(outformat=="svg"){
	# nam<-paste(fasta,"_synteny.svg",sep="")
	# svg(nam,width=11,  height = 1.1*le)
  # }
	# for( i in 1:le){
		# if(out[[i]][1,11]==-1){
		# temp<-out[[i]]
		
		# temp[,11]<-rep(1,nrow(temp))
		
		# ma<-max(temp[,8])
		# plus<-which(as.character(temp[,1])=="+")
		# minus<-which(as.character(temp[,1])=="-")
		# ord<-rep(NA,nrow(temp))
		# levels(temp[,1])<-c("+","-")
		# if(length(plus)>0){
			# temp[plus,1]<-"-"
			# ord[plus]<-"-"
		# }
		# if(length(minus)>0){
			# temp[minus,1]<-"+"
			# ord[minus]<-"+"
		# }
			# for(j in 1:nrow(temp)){
				# temp[j,8]<-ma-out[[i]][j,7]
				# temp[j,7]<-ma-out[[i]][j,8]
				# temp[j,10]<-ma-out[[i]][j,9]
				# temp[j,9]<-ma-out[[i]][j,10]
				# }
			# #out[[i]]<-temp	
		# }
		
	# }
	
	# #mas<-out[[1]][1,9]
	# mas<-wind
	 # for(i in 1:le){
		 # te<-out[[i]]
		 # s2<-mas-as.numeric(te[1,9])
		# for(ii in 1:nrow(te)){
        
         
		 # #print(as.numeric(te[1,9])[1])
        # #s<-c(s,s2+mas)
         # #se<-c(se,s2+mi)
         # te[ii,7]<-as.numeric(te[ii,7])+s2
         # te[ii,8]<-as.numeric(te[ii,8])+s2
         # te[ii,9]<-as.numeric(te[ii,9])+s2
         # te[ii,10]<-as.numeric(te[ii,10])+s2
         
		# }
		# out[[i]]<-te
     # }
	
  # #
  
  
  
  # d_vect<-c()
  # s_vect<-c()
  
	# plot(1,1, type="n",xlim=c(1,wind*2),ylim=c(1,count),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
    # for(i in 1:le){
    
        # temp<-out[[i]]
		# nam3<-names(out)[i]
        # mi<-min(as.numeric(temp[,7]))
        # ma<-max(as.numeric(temp[,8]))

        # lines(c(mi,ma),c(1+count-i*4,1+count-i*4))
        # lines(c(mi,ma),c(2+count-i*4,2+count-i*4))
        # for(j in 1:nrow(temp)){

          # m<-(as.numeric(temp[j,8])-as.numeric(temp[j,7]))/2+as.numeric(temp[j,7])
          # mm<-as.numeric(temp[j,8])-as.numeric(temp[j,7])
          # n<--1
          # if(temp[j,1]=="+"){
            # n<-1
          # }
          # if(n==-1){
            # n<-0
          # }
          # color<-"white"
          # tcolor<-na.omit(grep(temp[j,5],cdhit_result))
          # tcolor2<-na.omit(grep(temp[j,4],cdhit_result))
          # tcolor<-c(tcolor,tcolor2)
          # if(length(tcolor)>0){
             # color<-names(cdhit_result)[tcolor]
           # }
          # rect(as.numeric(temp[j,7]),0.5+n+count-i*4,as.numeric(temp[j,8]),1.5+n+count-i*4, col=color)
          # text(m,n+count-i*4+1+0.25,gsub("\"","",gsub("NA ","",gsub("orf.*","",temp[j,4]))),cex=0.4,font=4)
		  # text(m,n+count-i*4+1-0.25,gsub("\"","",temp[j,5]),cex=0.4,font=4)
        # }
        # n<-as.numeric(temp[1,11])
        # if(n==-1){
          # n<-0
        # }
        # rect(as.numeric(temp[j,9]),0.5+n+count-i*4,as.numeric(temp[j,10]),1.5+n+count-i*4, col=2)
        # text(wind,count-i*4,paste(nam3, temp[1,12], sep=" "), cex=0.8, font=4)
  	
    # }
 # dev.off()	

# }

#plot_function4(out2[[1]],cd,wind=synteny_window,outformat="pdf", fasta=name)

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

# # function to exclude very similar organism based on a phylogentic tree to reduce complexity
# exclude_similars<-function(dis, thres=0.01,ooi){
# i<-1

# o<-grep(ooi, rownames(dis))
# nam<-colnames(dis)[o]
	# temp<-which(dis[,o]<=thres)
	# nam<-na.omit(match(nam, rownames(dis)))
	# nam<-na.omit(match(nam, temp))
	# if(length(nam)>0){
		# temp<-temp[-nam]
	# }
	
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		
	# }
	



# while(nrow(dis)>i){
	
	# nam<-colnames(dis)[i]
	# temp<-which(dis[,i]<=thres)
	# nam<-na.omit(match(nam, rownames(dis)))
	# nam<-na.omit(match(nam, temp))
	# if(length(nam)>0){
		# temp<-temp[-nam]
	# }
	
	# if(length(temp)>0){
		# dis<-dis[-temp,-temp]
		
	# }
	# i<-i+1
# }
# dis
# }

exclude_similars<-function(dis, thres=0.01,ooi){
i<-1


for(j in 1:length(ooi)){
	o<-match(ooi, gsub("\\..*","",rownames(dis)))
	nam_ref<-colnames(dis)[o]
	temp<-which(dis[,o[j]]<=thres)
	
	temp<-setdiff(temp,o)
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		
	}
	}



while(nrow(dis)>i){
	o<-match(ooi, gsub("\\..*","",rownames(dis)))
	nam<-colnames(dis)[i]
	o2<-match(nam,colnames(dis))
	temp<-which(dis[,i]<=thres)
	temp<-setdiff(temp,c(o2,o))
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		
	}
	i<-i+1
}
dis
}



# Execute functions



fasta<-readLines(filename)
#fasta<-as.character(fasta[,1])
#fasta<-sub("/[0123456789]*.[0123456789]*","",fasta)
coor<-export_ncRNA_coordinates(fasta)
coor<-remove_overlapping_homologs(coor)

#seq1<-coor[[2]]
#coor<-coor[[1]]

#coor[,"sequence"]<-gsub("-","",coor[,"sequence"])

# fasta3<-gsub("-","",seq1[2])
# fasta2<-c(seq1[1],fasta3)
# for(j in 1:nrow(coor)){
	# fasta2<-c(fasta2,as.character(coor[j,"Full_header"] ))
	# fasta2<-c(fasta2, as.character(coor[j,"sequence"]))
# }
# write.table(fasta2, file=filename, row.names=F, col.names=F, quote=F)


#coor2<-coor
#taxi<-gsub(".*taxID:","",coor[,"Full_header"])

# assign matching Refseq ID based on the Accesion Number
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
na<-which((coor[,"fin"]=="NA"))
if(length(na)>0){
	coor<-coor[-na,]
}


#if more than 1 homolog is detected for one organism or the same Refseq ID, keep only the homolog with the highest identity to the input
ooi_pos<-grep(ooi,coor[,"fin"])

tempfasta<-tempfile()
tempfasta2<-tempfile()
tmpfas<-c()
for(i in 1:nrow(coor)){
	tmpfas<-c(tmpfas,coor[i,"Full_header"],coor[i,"sequence"])
}
writeLines(tmpfas,con=tempfasta)
mafft(filename=tempfasta,outname=tempfasta2,mode="fast")
dat<-read.phyDat(tempfasta2, format="fasta", type="DNA")

dm <- dist.hamming(dat)
dm<-as.matrix(dm)
iden<-1-dm[,ooi_pos]

coor<-cbind(coor,iden)


coor<-coor[order(iden, decreasing=T), ]




# if(mixed_sources==T){
	# mafft(filename=filename,outname="aligned_sRNA.fasta")
	# dat<-read.phyDat("aligned_sRNA.fasta", format="fasta", type="DNA")
	# unlink("aligned_sRNA.fasta")
	# dm <- dist.hamming(dat)
	# dm<-as.matrix(dm)
	# iden<-dm[,1]
	
	# coor<-cbind(coor,iden)
	# coor<-coor[order(iden, decreasing=T), ]
	
	# # remove sRNAs from the same organisms with overlapping coordinates
	# orgs<-unique(coor[,1])
	# over<-c()
	# for(i in orgs){
		# tmp<-which(coor[,1]==i)
		# if(length(tmp)>1){
			# coordina<-function(coor){
				# out<-as.numeric(coor[3]):as.numeric(coor[4])
				# out
			# }
			# coordi<-apply(coor[tmp,],1,coordina)
			# if(is.matrix(coordi)==T){
				# coordi<-as.list(data.frame(coordi))
			# }
			
			# j<-1
			# while(length(coordi)>1 & j<length(coordi)){
			
				# tmp<-c()
				# for(jj in (j+1):length(coordi)){
					# tmp<-c(tmp, length(intersect(coordi[[j]],coordi[[jj]])))
				# }
				# tmp<-which(tmp>0)
				
				# if(length(tmp)>0){
				
					# over<-c(over, names(coordi)[tmp])
					# coordi<-coordi[-tmp]
					
				# }
				# j<-j+1
			# }
			
			
			
			
		# }
	
	
	# }
	# if(length(over)>0){
		# over<-na.omit(match(over, rownames(coor)))
		# if(length(over)>0){
			# coor<-coor[-over,]
		# }
	# }
# }

# if(refseq_required==F){
	# if(is.matrix(coor)==T){
			# na<-which(is.na(coor[,"fin"]))
			# if(length(na)>0){
				# coor[na,"fin"]<-coor[na,"Accesion_number"]
			# }
		# }
# }

# # remove entries without Refseq ID
# if(refseq_required==T){
	if(is.matrix(coor)==T){
		na<-which(is.na(coor[,"fin"]))
		if(length(na)>0){
			coor<-coor[-na,]
		}
	}
#}


keep<-c(ooi,wildcard)
keep<-na.omit(match(keep,gsub("\\..*","",coor[,"fin"])))
keep2<-coor[keep,]
coor<-coor[-keep,]

# remove entries with duplicated Refseq ID and entries without matching Refseq ID
#if(duplicates_allowed==F){
	
		dup<-which(duplicated(coor[,"fin"]))
		if(length(dup)>0){
			coor<-coor[-dup,]
			#coor2<-coor2[-dup,]
		}
	
	
#}
coor<-rbind(keep2,coor)
# keep only homologs represented in the coprarna reference file	

#if(coprarna_compatible==T){
	copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
		se<-function(x){
			out<-grep(x, copref[,1])[1]
			if(length(out)==0){
				out<-NA
			}
			out
		}
	notinlist<-which(is.na(unlist(lapply(gsub("\\..*","",coor[,"fin"]),se))))
	if(length(notinlist)>0){
		coor<-coor[-notinlist,]	
	}
#}

# write fasta file with all sequences that have a assigned Refseq ID and are available for CopraRNA

# fasta3<-c()
		# for(i in 1:nrow(coor)){
			# fasta3<-c(fasta3, paste(">",coor[i,"fin"],"|",gsub(">","",coor[i,"Full_header"]),sep=""))
			# fasta3<-c(fasta3, as.character(coor[i,"sequence"]))
			
		# }
	

# write.table(fasta3, file=outfile, row.names=F, col.names=F, quote=F)


# create short names from organism name and Refseq ID

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:nrow(coor)){
	tnam<-grep(gsub("\\..*","",coor[i,"fin"]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}

nam2<-c()
for(i in 1:length(nam)){
	temp1<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp1<-paste(temp1,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp1<-paste(temp1, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp1)
}
nam2<-paste(nam2,coor[,"fin"], sep="_")

coor<-cbind(coor,nam2)
# extract information of sRNA gene neighborhood for synteny analysis and draw synteny pdf

#save(coor, file="coor.Rdata")

orgs<-unique(coor[,1])

command<-paste("python3 ", script_path, " -s ", db_path, " -rRNA ", paste(orgs, collapse=" "))
print(command)
rRNA<-system(command, intern=T)

startp<-grep(">",rRNA)[1]
rRNA<-rRNA[startp:length(rRNA)]

empty<-grep("No 16sRNA sequence found",rRNA)
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

	
tempfas<-tempfile()
tempfas2<-tempfile()
writeLines(rRNA,con=tempfas)
#mafft(filename=tempfas,outname=tempfas2)
command<-paste("mafft --thread 30 --retree 2 --maxiterate 0  --quiet --inputorder ", tempfas, " > ", tempfas2, sep="" )
system(command)


dat<-read.phyDat(tempfas2, format="fasta", type="DNA")
#dm <- dist.ml(dat, model="F81")
#dm2<-as.matrix(dm)
dm <- dist.hamming(dat)
dm<-as.matrix(dm)


ma<-match(rownames(dm), coor[,1])
rownames(dm)<-gsub("\\..*","",coor[ma,"fin"])
colnames(dm)<-rownames(dm)
dm_red<-exclude_similars(dm,ooi=c(ooi,wildcard),thres=0.02)
dim(dm_red)

get_similars<-function(dm, ref_orgs, thres=0.01, num=rep(3,length(ref_orgs)),cands){
	outp<-c()
	for(i in 1:length(ref_orgs)){
		pos<-match(ref_orgs[i], colnames(dm))
		
		tmp<-c()
		count<-1
		
		while(length(tmp)<num[i] & count<=100){
		tmp<-which(dm[,pos]<=thres)		
		tmp<-dm[pos,tmp]
		ex<-na.omit(match(names(tmp),cands))
			if(length(ex)>0){
				tmp<-tmp[-ex]
			}
		thres<-thres+0.002	
		count<-count+1
		}
			if(length(tmp)>=num[i]){
				le<-length(tmp)
				le2<-max(floor(le/num[i]),1)
				le3<-seq(le2,le,by=le2)
				outp<-c(outp,names(tmp)[le3])
			} else {
				outp<-c(outp,names(tmp))
			}
			
	}
	outp
}
# S. oneidensis, V. cholerae, Aeromonas hydrofila, Y. pseudotubecolosis YPIII, S. enterica, Xenorhabdus Nematophila, Aliivibrio salmonicida
cnds<-get_similars(dm, ref_orgs=c(ooi,wildcard), thres=0.015, num=c(4,4,4,7,4,3,4), cands=rownames(dm_red))
cnds<-na.omit(c(colnames(dm_red),cnds))
length(unique(cnds))

n2<-unique(match(cnds,gsub("\\..*","",coor[,"fin"])))

orgs<-(coor[n2,1])

orgs<-c(orgs,"CP053706.1")

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


tempfas<-tempfile()
tempfas2<-tempfile()
writeLines(rRNA,con=tempfas)
#mafft(filename=tempfas,outname=tempfas2)
command<-paste("mafft --thread 40 --maxiterate 1000 --localpair --quiet --inputorder ", tempfas, " > ", tempfas2, sep="" )
system(command)


dat2<-read.phyDat(tempfas2, format="fasta", type="DNA")
dm3 <- dist.ml(dat2, model="F81")
dm4<-as.matrix(dm3)


#save(dm2, file="distances.Rdata")
	
treeNJ <- NJ(dm3)

	#tree <- pratchet(dat)          # parsimony tree
	#tree <- nnls.phylo(tree, dm)

	#fitStart = pml(treeNJ, dat, k=4)
	#fitJC = optim.pml(fitStart, model="GTR", optGamma=TRUE, rearrangement="none")
	#treeNJ<-tree
#cd Cop
cop<-match(unique(cnds),gsub("\\..*","",coor[,"fin"]))

fast<-c()
for(i in cop){
	fast<-c(fast,paste0(">",gsub("\\..*","",coor[i,"fin"])),coor[i,"sequence"])
}
writeLines(fast,con="~/media/jens@margarita/Copra2_paper/Glassgo/RyhB/copra_inp.fasta")
