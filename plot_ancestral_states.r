require(ggplot2)
require(ggtree)
require(ggimage)
require(ape)


#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/plot_ancestral_states.r 

piecolors=c("#b0b0b0" ,"#799943", "#5680b0", "#9c6488" ,"#e4c054" ,"#e26a5f" ,"#b43768" ,"#006362" ,"#614736", "#294258", "#77b1b5")
load("ancestral_site_recons.Rdata")
aligned=T

pdf("ancestral_reconstruction.pdf",useDingbats=FALSE)
for(i in 1:length(anc_res)){
	tree<-anc_res[[i]][[3]]
	
	ab<-grep("CP014266",tree$tip.label)
	if(length(ab)>0){
		tree$tip.label[ab]<-"Aci_baumannii_CP014266"
	}
	inner<-unique(tree[[1]][,1])

	opts<-anc_res[[i]][[6]]
	subs<-anc_res[[i]][[7]]

	o1<-sort(unique(opts))
	s1<-sort(unique(subs))


	#pos<-na.omit(match(gsub(".*_","",names(opts)),gsub(".*_","",tree$tip.label)))
	################

	pos2<-na.omit(match(gsub(".*_","",tree$tip.label),gsub(".*_","",names(opts))))
	opts<-opts[pos2]
	pos<-na.omit(match(gsub(".*_","",names(opts)),gsub(".*_","",tree$tip.label)))

	d1<-matrix(,length(opts),length(o1))
	#colnames(d1)<-o1
	rownames(d1)<-names(opts)
	for(j in 1:length(o1)){
		temp<-which(opts==o1[j])
		d1[temp,j]<-"A"

	}
	#o1<-o1[which(o1>0)]
	pie1<-c("white","#b0b0b0" ,"#799943", "#5680b0", "#9c6488" ,"#e4c054" ,"#e26a5f" ,"#b43768" ,"#006362" ,"#614736", "#294258", "#77b1b5")
	pie1<-pie1[o1+1]

	d1<-data.frame(lab=tree$tip.label[pos], d1)

	#####################

	pos2<-na.omit(match(gsub(".*_","",tree$tip.label),gsub(".*_","",names(subs))))
	subs<-subs[pos2]
	pos<-na.omit(match(gsub(".*_","",names(opts)),gsub(".*_","",tree$tip.label)))

	d2<-matrix(,length(subs),length(s1))
	colnames(d2)<-s1
	rownames(d2)<-names(subs)
	for(j in 1:length(s1)){
		temp<-which(subs==s1[j])
		d2[temp,j]<-"A"

	}
	#s1<-s1[which(s1>0)]
	pie2<-c("white","#b0b0b0" ,"#799943", "#5680b0", "#9c6488" ,"#e4c054" ,"#e26a5f" ,"#b43768" ,"#006362" ,"#614736", "#294258", "#77b1b5")
	pie2<-pie2[s1+1]
	colnames(d2)<-paste(1:ncol(d2),"sub",sep="_")

	d2<-data.frame(lab=tree$tip.label[pos], d2)

	
	
	
	dat<-data.frame(anc_res[[i]][[1]]$states, node=sort(inner))
	sites<-as.numeric(colnames(anc_res[[i]][[1]]$states))
	piecolors2<-piecolors[sites]
	if(aligned==T){
	n1<-ggtree(tree) + geom_tiplab(size=2,hjust=-0.04, align = TRUE, linetype = 'dotted', linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6)) + ggtitle (names(anc_res)[i])
	a<-ggplot_build(n1)
	dis_a<-a[[1]][[3]][1,1]
	dis_b<-a[[1]][[4]][1,1]
	dis_c<-dis_b-dis_a
	command1<-paste("p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab(hjust=-",20*dis_c,", size=2,align = TRUE, linetype = 'dotted', linesize = 0.05) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))",sep="")

	} else {
		n1<-ggtree(tree) + geom_tiplab(size=2,hjust=-0.04, linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))
		a<-ggplot_build(n1)
		dis_a<-a[[1]][[2]][1,2]
		dis_b<-a[[1]][[3]][1,2]
		dis_c<-dis_b-dis_a
		command1<-paste("p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab(hjust=-",15*dis_c,", size=2, linesize = 0.05) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))",sep="")

	}
	
	#command1<-paste("p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab(hjust=-",10*dis_c,", size=2,align = TRUE, linetype = 'dotted', linesize = 0.05) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))",sep="")
	#command1<-paste("p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab( size=2,align = TRUE, linetype = 'dotted', linesize = 0.05) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))",sep="")
	#command1<-"p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab(size=2,hjust=-",dis_a*36,, align = TRUE, linetype = 'dotted', linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))"
	for(j in 1:(ncol(d1)-1)){
			tmp<-paste(" + geom_point(aes(shape=X",j,",x=",dis_a,"),size=3,na.rm=T,color=pie1[",j,"],show.legend = F)",sep="")
			command1<-paste(command1,tmp,sep="")
	}
		
	for(j in 1:(ncol(d2)-1)){
			tmp<-paste(" +  geom_point(aes(shape=X",j,"_sub,x=",dis_a+5*dis_c,"),size=3,na.rm=T,color=pie2[",j,"],show.legend = F)",sep="")
			command1<-paste(command1,tmp,sep="")
	}	
	
	
	#pdf(paste(names(anc_res)[i],"_node_ref.pdf",sep=""),useDingbats=FALSE)
	plot(tree,cex=0.4, main =names(anc_res)[i])
	nodelabels(cex=0.5, frame="none")
	#dev.off()
	#pdf(paste(names(anc_res)[i],"_ancestral_states.pdf",sep=""),useDingbats=FALSE)
	eval(parse(text=command1))
	#p<-ggtree(tree) + geom_tiplab(size=1.5,hjust=-0.05, align = TRUE, linetype = "dotted", linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6)) 

	pies <- nodepie(dat, cols=1:(ncol(dat)-1), color=piecolors2[1:(ncol(dat)-1)],alpha=1)
	p1<-inset(p, pies, width=.1, height=.1)
	plot(p1)
	#dev.off()
	

}

dev.off()






# dat<-data.frame(anc_res[[i]][[1]]$states, node=sort(inner))

# n1<-ggtree(tree) + geom_tiplab(size=1.5,hjust=-0.05, align = TRUE, linetype = 'dotted', linesize = 0.05,)
# a<-ggplot_build(n1)

# dis_a<-a[[1]][[3]][1,1]
	

# command1<-"p<-ggtree(tree) %<+% d1 %<+% d2 + geom_tiplab(size=2,hjust=-0.4, align = TRUE, linetype = 'dotted', linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6))"
# #command1<-"p<-ggtree(tree) %<+% d1 + geom_tiplab(size=1.5,hjust=-0.05, linetype = 'dotted', linesize = 0.05,) "

# for(j in 1:(ncol(d1)-1)){
		# #tmp<-paste(" + geom_point(aes(shape=X",j,",x=x+",dis_a*5,"),size=3,na.rm=T,color=pie1[",j,"],show.legend = F)",sep="")
		# tmp<-paste(" + geom_point(aes(shape=X",j,",x=",dis_a+0.05,"),size=5,na.rm=T,color=pie1[",j,"],show.legend = F)",sep="")
		# command1<-paste(command1,tmp,sep="")
	# }
	
# for(j in 1:(ncol(d2)-1)){
		# #tmp<-paste(" +  geom_point(aes(shape=X",j,"_sub,x=x+",dis_a*10,"),size=3,na.rm=T,color=pie2[",j,"],show.legend = F)",sep="")
		# tmp<-paste(" +  geom_point(aes(shape=X",j,"_sub,x=",dis_a+0.12,"),size=5,na.rm=T,color=pie2[",j,"],show.legend = F)",sep="")
		# command1<-paste(command1,tmp,sep="")
	# }	
# eval(parse(text=command1))
# pdf("polA_anc.pdf", useDingbats=F)
# #p<-ggtree(tree) + geom_tiplab(size=1.5,hjust=-0.05, align = TRUE, linetype = "dotted", linesize = 0.05,) +  coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6)) 

# pies <- nodepie(dat, cols=1:(ncol(dat)-1), color=piecolors[1:(ncol(dat)-1)],alpha=1)
# p1<-inset(p, pies, width=.11, height=.11)
# p1
# dev.off()

