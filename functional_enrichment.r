require(genbankr)
require(topGO)


#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/functional_enrichment.r --args num=200

# get absolute path  
#path="/home/jens/CopraRNA-git/coprarna_aux/"
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("functional_enrichment.r","",path)
print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")

load("copra_results_all.Rdata")
all_orgs<-names(copra_results)

num=200   # number of top CopraRNA predictions used for the enrichment


# transforming arguments into valid variables 
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

enrich_list<-vector("list",length(all_orgs))
names(enrich_list)<-all_orgs


for(i in 1:length(all_orgs)){
print(i)
cop_all<-as.matrix(copra_results[[i]])
present<-which(cop_all[,all_orgs[i]]!="")
cop_all<-cop_all[present,]
e<-grep("Annotation", colnames(cop_all))

cop<-cop_all[,4:(e-1)]

# int<-cop[,all_orgs[i]]
# int<-strsplit(int, "\\|")
# int<-as.numeric(do.call(rbind,int)[,3])
# ab<-which(int>=int_thres)
# if(length(ab)>0){
	# cop<-cop[-ab,]
# }
gbk<-paste("~/For_CopraRNA2.0/Genomes/",all_orgs[i],".gb.gz",sep="")
gbk1<-readLines(gbk)
org<-grep("ORGANISM",gbk1)
org<-gbk1[org]
org<-gsub(".*ORGANISM. ","",org)
org<-gsub(" .*","",org)
gb<-readGenBank(gbk)

cd<-cds(gb)
cd<-GenomicRanges::mcols(cd)

url <- "https://www.uniprot.org/uploadlists/"

cop_ids<-unname(cop[,all_orgs[i]])
cop_ids<-gsub("\\(.*","",cop_ids)


id_pos<-na.omit(match(tolower(cop_ids), tolower(cd$locus_tag)))

ids<-unlist(cd$locus_tag)[id_pos]
ids_keep<-ids


id_map<-cbind(ids,ids)
if(all_orgs[i]=="NC_002505"){
	ids<-gsub("VC","VC_",ids)
	id_map[,1]<-ids
}
query_string=paste(na.omit(ids), collapse=" ")
#get go_terms, uniprot id, go_ids
params = list(
  from = "GENENAME",
  to = "ACC",
  columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
  format = "tab",
  query = query_string
)

r <- httr::POST(url, body = params, encode = "form")
anno<-httr::content(r,encoding ="UTF-8")
if(is.null(anno)){
	ids<-unlist(cd$old_locus_tag)[id_pos]	
	id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
	query_string=paste(na.omit(ids), collapse=" ")
	#get go_terms, uniprot id, go_ids
	params = list(
	  from = "GENENAME",
	  to = "ACC",
	 
	  columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
	  format = "tab",
	  query = query_string
	)

	r <- httr::POST(url, body = params, encode = "form")
	anno<-httr::content(r,encoding ="UTF-8")
}


if(is.null(anno)){
	j<-1
	r2<-"\n"
	while(r2 == "\n" & j <= 10){
		r <- httr::GET(paste("http://rest.kegg.jp/find/genes/",na.omit(unlist(cd$locus_tag))[j] ,sep=""))
		r2<-httr::content(r,encoding ="UTF-8")
		j<-j+1
		if(is.null(r2)){
			r2<-"\n"
		}
	}
	
	if(r2 == "\n"){
		
		j<-1
		while(r2 == "\n" & j <= 10){
		#print(r2)
			r <- httr::GET(paste("http://rest.kegg.jp/find/genes/",na.omit(unlist(cd$old_locus_tag))[j] ,sep=""))
			r2<-httr::content(r,encoding ="UTF-8")
			j<-j+1
			if(is.null(r2)){
				r2<-"\n"
			}
		}	
		if(r2 != "\n"){
			r2<-gsub(":.*","",r2)
			ids<-unlist(cd$old_locus_tag)[id_pos]	
			ids<-paste(r2,":",ids,sep="")
			id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
			query_string=paste(na.omit(ids), collapse=" ")
			#get go_terms, uniprot id, go_ids
			params = list(
			  from = "KEGG_ID",
			  to = "ACC",
			   columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
			  format = "tab",
			  query = query_string
			)

			r <- httr::POST(url, body = params, encode = "form")
			anno<-httr::content(r,encoding ="UTF-8")
		} else {
			next
		}
	} else {
		r2<-gsub(":.*","",r2)
		ids<-unlist(cd$locus_tag)[id_pos]	
		ids<-paste(r2,":",ids,sep="")
		id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
		query_string=paste(na.omit(ids), collapse=" ")
		#get go_terms, uniprot id, go_ids
		params = list(
		  from = "KEGG_ID",
		  to = "ACC",
		   columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
		  format = "tab",
		  query = query_string
		)

		r <- httr::POST(url, body = params, encode = "form")
		anno<-httr::content(r,encoding ="UTF-8")
	}
}

if(is.null(anno)==F ){
write.table(anno,file="tmp_uniprot.txt",row.names=F, col.names=F, quote=F)
anno<-read.csv("tmp_uniprot.txt", header=F, skip=1, sep="\t")
if(is.null(anno)==F & nrow(anno) > 100){
path<-gsub(",","",anno[,7])
path<-gsub(" ","_",path)
path<-gsub("PATHWAY:_","",path)

term<-anno[,2]
term<-gsub(" ","",term)
term<-gsub(";;","",term)
term<-gsub(";",", ",term)

mat<-match(as.character(anno[,ncol(anno)]), id_map[,1])
an<-id_map[mat,2]

coprarank<-match(toupper(an),gsub("\\(.*","",toupper(unlist(cop[,all_orgs[i]]))))
term<-cbind(toupper(an),term)

term2<-cbind(anno,id_map[mat,],coprarank)
write.table(term, file="enrich_anno.txt", sep="\t", quote=F, row.names=F,col.names=F)
geneID2GO <- readMappings("enrich_anno.txt")
geneNames <- names(geneID2GO)
#myInterestingGenes<-cop_ids[1:num]
#geneList<-as.numeric(cop_all[,2])
geneList<-1:nrow(cop_all)
names(geneList)<-toupper(cop_ids)


geneList<-na.omit(geneList)
#geneList<-na.omit(geneList)

#geneNames <- names(geneID2GO)
#myInterestingGenes <- geneNames[1:200]
#geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
#names(geneList) <- geneNames


topDiffGenes <- function(allScore) { 
	#return(allScore < 0.01) 
	return(allScore < (num+1))
}
genesOfInterest<-topDiffGenes(geneList)

# topDiffGenes2 <- function(allScore) { 
	# return(allScore < 0.2) 
	# #return(allScore < 201)
# }
# genesOfInterest2<-topDiffGenes2(geneList)


#geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
#names(geneList) <- geneNames
go_cat<-c("BP","MF","CC")
temp_enrich<-list()
temp_go_dat<-list()
nodes<-c()
for(jj in go_cat){
GOdata <- new("topGOdata",ontology = jj,allGenes = geneList,geneSel = topDiffGenes,annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 2)
temp_go_dat[[jj]]<-GOdata
#result01FS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
#result01ks <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
classicFS <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#classicKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
#nodes<-c(nodes,result01ks@geneData[4])
#nodes<-c(nodes,result01FS@geneData[4])
nodes<-c(nodes,classicFS@geneData[4])
#nodes<-c(nodes,classicKS@geneData[4])
#GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)


allRes <- GenTable(GOdata,classicFS=classicFS ,ranksOf = "classicFS", topNodes = min(nodes)) # classicKS = classicKS , weight01FS=result01FS, weight01ks=result01ks


allRes$genes <- sapply(allRes$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata, x) 
  paste(genes[[1]][genes[[1]] %in% names(which(genesOfInterest))],collapse=",")
})

allRes$ranks <- sapply(allRes$genes, function(x)
{
	if(length(x)>0){
		x<-strsplit(x,",")[[1]]
		posi<-paste(match(x, toupper(gsub("\\(.*","",cop[,all_orgs[i]]))), collapse=",")
		posi
  }
})


#allRes <- GenTable(GOdata,   ks=resultKs,  topNodes = 1000)
#write.table(allRes, file="enrichment.txt", sep="\t")
temp_enrich[[jj]]<-allRes
}

enrich_list[[i]][[1]]<-temp_enrich
enrich_list[[i]][[2]]<-temp_go_dat
enrich_list[[i]][[3]]<-term2
enrich_list[[i]][[4]]<-anno
enrich_list[[i]][[5]]<-cd
}
}
}
save (enrich_list, file="enrichment.Rdata")


