pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools", "metagenomeSeq","minpack.lm", "Hmisc", "reltools", "stats4","venn")

lapply( pkg, library, character.only = TRUE )
set.seed( 131 )
rm(list = ls()[grep("*.*", ls())])
#source("/home/hassani/scripts/sncm.fit_function.r")
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")
Palette.shape.1<-c(21, 21, 16)
Palette.shape.2<-c(13,17,19,15,25,10)
alpha=0.05
minZ=1.5
maxZ=8.5

DT<-list()
disp<-list()
mci<-list()
gp<-list()
theme_change <- theme(
#	legend.position="none",
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.border = element_rect(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.background = element_rect( fill="white", colour ="black" )
	)

	Bee<-c( 'Bee1','Bee2','Bee3','Bee4','Bee5','Bee6','Bee7','Bee8' )

	T1D1<-c('T1D11','T1D12','T1D13','T1D14','T1D15','T1D16')
	T1D3<-c('T1D31','T1D32','T1D33','T1D34','T1D35','T1D36')
	T1D5<-c('T1D51','T1D52','T1D53','T1D54','T1D55','T1D56')
	
	T2D1<-c('T2D11','T2D12','T2D13','T2D15','T2D16')
	T2D3<-c('T2D31','T2D32','T2D33','T2D34','T2D35','T2D36')
	T2D5<-c('T2D51','T2D52','T2D53','T2D54','T2D55','T2D56')
	
	T1<-c( T1D1,T1D3,T1D5 )
	T2<-c( T2D1,T2D3,T2D5 )
	
	All<-c( Bee,T1,T2 )

	All.S<-list(Bee,T1D1,T1D3,T1D5,T2D1,T2D3,T2D5)

# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T )
	mat=as.matrix( mat )
	OTU=otu_table( mat, taxa_are_rows=T ) 
	tax=read.table( "taxonomy.txt", sep="\t", row.names=1, header=1 )
	tax=as.matrix( tax )
	TAXA=tax_table( tax )
	sd=read.table( "sample_data.txt", sep="\t", row.names=1, header=1 )
	SD=sample_data( sd )

	physeq=phyloseq(OTU, TAXA, SD) 
	physeq.MC=subset_taxa( physeq, !Genus %in% c("Mitochondria","Chloroplast") )
	physeq_B=prune_samples(sample_sums(physeq.MC)>=1000, physeq.MC)

	tax.B=tax_table(physeq_B)
	Sd=sample_data(physeq_B)

# transformt counts to rel. abundance
	physeq.B.ra = transform_sample_counts(physeq_B, function(x) x/sum(x))

# computing arithmetic mean of rel. abund. by genotype
	DT.B.ra=data.table(otu_table(physeq.B.ra), keep.rownames=T, key="rn")
	DT.B.taxa <- data.table(as.data.frame(tax.B), keep.rownames=T, key="rn")
	ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
	DT.B.taxa$ColPhy <- ifelse(!DT.B.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.B.taxa$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT.B.taxa$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.B.taxa$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))
	DT.B.m=merge(DT.B.ra, DT.B.taxa)
	OTU_B=t(otu_table(physeq_B))
	TAXA_B=data.frame(DT.B.taxa)	
	rownames(TAXA_B)<-TAXA_B$rn
	TAXA_B<-TAXA_B[,-1]

for (i in All.S) {
	DT.B<-DT.B.m[,c( "rn",..i )]
	DT.B$mean<-DT.B[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
	list.B=DT.B[DT.B$mean!=0]$rn
		
	print(paste0("Computing dispersal for ",i[[1]], " samples ..."))	

	spp=OTU_B[i,list.B]
	
	N <- mean(apply(spp, 1, sum))
	p.m <- apply(spp, 2, mean)
	p.m <- p.m[p.m != 0]
	p <- p.m/N
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]
	d = 1/N

	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001))
  	m.ci <- confint(m.fit, 'm', level=0.95)
	
	disp[[i[[1]]]]<-as.numeric(coef(m.fit))
   	mci[[i[[1]]]]=as.numeric(coef(m.fit)-m.ci[1])
	print(paste0("disparsal rate = ",coef(m.fit)))
	
	DT.C=data.table(C, key="Row.names")
	DT.C$x_y=DT.C$y*DT.C$x
	setnames(DT.C,"Row.names","rn")
	DT1=merge(DT.C, DT.B.taxa, by.x="rn", by.y="rn")
	DT.mean=DT.B[,c('rn','mean')]
	DT2=merge(DT1, DT.mean, by.x="rn",by.y="rn")
	DT2$names1=ifelse(DT2$y > .8, DT2$Genus,"")
	DT2$names2=ifelse(DT2$y > .8, DT2$ASV,"")
	
		
	pp1=ggplot( DT2, aes( x=mean,y=y ) )
	p1=pp1+geom_jitter( aes( color=ColPhy ), size=1 )+
	scale_color_manual( values=Palette.phylum ) +
	geom_hline(yintercept = .8, linetype="longdash")+
	scale_y_continuous( limits=c( 0,1.1 ),breaks=c( 0,0.2,0.4,0.6,0.8,1 ) ) +
	scale_x_continuous( limits=c( 0,0.55 ),breaks=c( 0,0.1,0.2,0.3,0.4,0.5 ) ) +
	theme_bw( )+theme_change

#	print( p1 )
	
	DT[[i[[1]]]]<-DT2[,c(1,3)]
	gp[[i[[1]]]]<-p1
	
#	pdf("dispersal_bacteria_all.pdf" ,useDingbats=FALSE)
#	gridExtra::grid.arrange(p1, nrow=2, ncol=2)
#	dev.off()
	
}	
print("##############################")
print("##############################")
print("##############################")

	
#	pdf("dispersal_bacteria.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(gp$T1D11,gp$T1D31,gp$T1D51, 
				gp$T2D11,gp$T2D31,gp$T2D51,
				gp$Bee1,nrow=3, ncol=3)	
#	dev.off()

	DT.Bee<-DT$Bee
	setnames(DT.Bee, "y","Bee")
	
	DT.T1D1<-DT$T1D11
	setnames(DT.T1D1, "y","T1D1")
	
	DT.T1D3<-DT$T1D31
	setnames(DT.T1D3, "y","T1D3")
	
	DT.T1D5<-DT$T1D51
	setnames(DT.T1D5, "y","T1D5")
	
	DT.T2D1<-DT$T2D11
	setnames(DT.T2D1, "y","T2D1")
	
	DT.T2D3<-DT$T2D31
	setnames(DT.T2D3, "y","T2D3")
	
	DT.T2D5<-DT$T2D51
	setnames(DT.T2D5, "y","T2D5")
	
	DT.T1=DT.T1D1[DT.T1D3][DT.T1D5]
	DT.T2=DT.T2D1[DT.T2D3][DT.T2D5]
			
	DT.all1=merge(DT.Bee, DT.T1, all=TRUE)
	DT.all=merge(DT.all1, DT.T2, all=TRUE)
	DT.all[is.na(DT.all)]<-0
	DT.all=merge(DT.all, DT.B.taxa)
	
	ASV.Bee=DT.all[DT.all$Bee>0.8]$ASV
	ASV.T1D1=DT.all[DT.all$T1D1>0.8]$ASV
	ASV.T1D3=DT.all[DT.all$T1D3>0.8]$ASV
	ASV.T1D5=DT.all[DT.all$T1D5>0.8]$ASV
	ASV.T2D1=DT.all[DT.all$T2D1>0.8]$ASV
	ASV.T2D3=DT.all[DT.all$T2D3>0.8]$ASV
	ASV.T2D5=DT.all[DT.all$T2D5>0.8]$ASV
	
	graph.T1=venn(list(ASV.T1D1, ASV.T1D3, ASV.T1D5), zcolor="style", snames=c("D1","D3","D5"))
	attr.T1<-attr(graph.T1,"intersection")
	
	graph.T2=venn(list(ASV.T2D1, ASV.T2D3, ASV.T2D5), zcolor="style", snames=c("D1","D3","D5"))
	attr.T2<-attr(graph.T2,"intersection")
	
	ASV.T1=unique( c( ASV.T1D1,ASV.T1D3,ASV.T1D5 ) )
	ASV.T2=unique( c( ASV.T2D1,ASV.T2D3,ASV.T2D5 ) )
	graph.T12B=venn(list(ASV.T1, ASV.T2, ASV.Bee), zcolor="style", snames=c("T1","T2","Bee"))
	attr.T12B<-attr(graph.T12B,"intersection")
	
		

#for (i in All) {
	DT.B<-DT.B.m[,c("rn",..All)]
	DT.B$mean<-DT.B[,.(rowMeans(.SD,na.rm = T)),.SDcols = All]
	list.B=DT.B[DT.B$mean!=0]$rn
		
	print(paste0("Computing dispersal for all samples ..."))	

	spp=OTU_B[All,list.B]
	
	N <- mean(apply(spp, 1, sum))
	p.m <- apply(spp, 2, mean)
	p.m <- p.m[p.m != 0]
	p <- p.m/N
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]
	d = 1/N

	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001))
  	m.ci <- confint(m.fit, 'm', level=0.95)
	
	disp[[i[[1]]]]<-as.numeric(coef(m.fit))
   	mci[[i[[1]]]]=as.numeric(coef(m.fit)-m.ci[1])
	print(paste0("disparsal rate = ",coef(m.fit)))
	
	DT.C=data.table(C, key="Row.names")
	DT.C$x_y=DT.C$y*DT.C$x
	setnames(DT.C,"Row.names","rn")
	DT1=merge(DT.C, DT.B.taxa, by.x="rn", by.y="rn")
	DT.mean=DT.B[,c('rn','mean')]
	DT2=merge(DT1, DT.mean, by.x="rn",by.y="rn")
	DT2$names1=ifelse(DT2$y > .8, DT2$Genus,"")
	
	
	list1=c(attr.T12B$'Bee')
	list2=c(attr.T12B$'T1:Bee')
	list3=c(attr.T12B$'T1:T2:Bee')
	list4=c(attr.T12B$'T1:T2')
	list6=c(attr.T12B$'T1')
	
	list5=c( list1,list2,list3,list4,list6 )
	
	
	DT2$shape=ifelse( DT2$ASV %in% list1,"Bee",ifelse( DT2$ASV %in% list2, "T1B", ifelse( DT2$ASV %in% list3,"core", ifelse( DT2$ASV %in% list4,"T12", ifelse(DT2$ASV %in% list6, "T1", "else" )))))
	
	DT2$names2=ifelse(DT2$ASV %in% list5 | DT2$y > .8 , DT2$Genus,"")
		
	pp2=ggplot( DT2, aes( x=log2(mean),y=y ) )
	p2=pp2+geom_jitter( aes( color=ColPhy, shape=shape ), size=3 )+
	geom_text( aes( label=names2 ), size=2 )+
	scale_color_manual( values=Palette.phylum ) +
	scale_shape_manual( values=Palette.shape.2 )+
	geom_hline(yintercept = .8, linetype="longdash")+
	#scale_y_continuous( limits=c( 0,1.1 ),breaks=c( 0,0.2,0.4,0.6,0.8,1 ) ) +
#	scale_x_continuous( limits=c( 0,0.55 ),breaks=c( 0,0.1,0.2,0.3,0.4,0.5 ) ) +
	theme_bw( )+theme_change

#	pdf("frq_relabun.pdf" ,useDingbats=FALSE)
	print( p2 )
#	dev.off()
#}	
print("##############################")
print("##############################")
print("##############################")

		






















