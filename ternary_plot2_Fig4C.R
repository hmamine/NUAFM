pkg=c( "plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools","ggtern","PMCMRplus","venn", "metagenomeSeq","UpSetR", "venn","seqtime" )

lapply( pkg, library, character.only = TRUE )

rm( list = ls( )[ grep( "*.*", ls( ) ) ] )
set.seed( 131 )

#Palette.phylum <-c( "#EF5656","#88888a","#2BB065" )
#Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")
Palette.phylum <-c("#EF5656","#F7A415","#88888a","#2BB065")
#Palette.phylum <-c( "#88888a","#2BB065" )

ColPhylum<-c( "Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria" )
Shape.sign<-c( 21,16,15 )
lines <-data.frame( x=c( 0.5,0,0.5 ),y=c( 0.5,0.5,0 ),z=c( 0,0.5,0.5 ), xend=c( 1,1,1 )/3, yend=c( 1,1,1 )/3, zend=c( 1,1,1 )/3 )


add_seg<-geom_segment( data=lines, aes( x,y,z, xend=xend, yend=yend, zend=zend ),color="grey", size=1,linetype="dotdash" )

theme_change <- theme(
#	legend.position="none",
	plot.background = element_blank( ),
	panel.grid.minor = element_blank( ),
	panel.grid.major = element_blank( ),
	)
DT<-vector( "list" )
DT.out<-list( )
RN<-list( )

#log2 fold change
	folch=1.5
#adjusted p-value
	alpha=0.05


# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T )
	mat=as.matrix( mat )
	OTU=otu_table( mat, taxa_are_rows=T ) 
	tax=read.table( "taxonomy.txt", sep="\t", row.names=1, header=1 )
	tax=as.matrix( tax )
	TAXA=tax_table( tax )
	sd=read.table( "sample_data.txt", sep="\t", row.names=1, header=1 )
	SD=sample_data( sd )
	TREE=read.tree( "tree.nwk" )

	physeq=phyloseq( OTU, TAXA, SD ) 
	
	physeq.MC=subset_taxa( physeq, !Genus %in% c( "Mitochondria","Chloroplast" ) )
	
	physeq.prune <- prune_samples( sample_sums( physeq.MC )>=1000, physeq.MC )
	physeq.subset <- subset_samples( physeq.prune, Cond != "T3" )
	
	DT.taxa=data.table( as( tax_table( physeq.subset ),"matrix" ), keep.rownames=T, key="rn" )
	DT.taxa$ColPhylum <- ifelse( !DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Actinobacteriota", "Actinobacteria", ifelse(DT.taxa$Phylum == "Bacteroidota","Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes","Firmicutes", "Proteobacteria" ) ) ) )
	
	
#trimming OTUs with low occurance
	otu.I=otu_table( physeq.subset )
	taxa.I=tax_table( physeq.subset )
	sd.I=sample_data( physeq.subset )
	flt=filterTaxonMatrix( otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE )
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff( 1:nrow( taxa.I ),flt$filtered.indices ),] 
	dummyTaxonomy=c( "k__","p__","c__","o__","f__","g__","s__" )
	TAXA.filt=rbind( taxa.I.filtered, dummyTaxonomy )
	rownames( TAXA.filt )[nrow( TAXA.filt )]="SUM"
	rownames( otu.I.filtered )[nrow( otu.I.filtered )]="SUM"
	NewOTU=otu_table( otu.I.filtered, taxa_are_rows = TRUE )
	NewTaxa=tax_table( TAXA.filt )

	physeq.trim=phyloseq( NewOTU, NewTaxa, sd.I )
	
#	physeq.I=subset_samples( physeq.trim, Day %in% c( "D3", "Bee" ) )
	physeq.I=physeq.trim
#computing DE ASVs using fitZig	
	physeq.i=subset_samples( physeq.I, Cond %in% c( "T1","T2" ) ) 
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cond
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cond, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE )	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CondT1	-	CondT2, levels = finalMod1 )
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable( fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames( DT.zig )
	DT.zig<-data.table( DT.zig, key="rn" )

	DT12=DT.taxa[ DT.zig, ]

	DT12$Col1 <- ifelse( DT12$adj.P.Val < alpha & abs(DT12$logFC) > folch, DT12$ColPhylum,"NS" )
	LIST.12=DT12[ DT12$Shp1 != "NS" ]$rn
#	DT12$Shp1 <- ifelse( DT12$adj.P.Val < alpha & DT12$logFC > folch, "T1", ifelse( DT12$adj.P.Val < alpha & DT12$logFC < -folch, "T2","NS" ) )
	

	
#computing DE ASVs using fitZig
	physeq.i=subset_samples( physeq.I, Cond %in% c( "T1","Bee" ) )	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cond
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cond, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE )	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CondT1	-	CondBee, levels = finalMod1 )
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable( fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames( DT.zig )
	DT.zig<-data.table( DT.zig, key="rn" )

	DT1B=DT.taxa[ DT.zig, ]

	DT1B$Col <- ifelse ( DT1B$adj.P.Val < alpha & abs(DT1B$logFC) > folch,  DT1B$ColPhylum,"NS" )
	LIST.1B=DT1B[ DT1B$Col != "NS" ]$rn
#	DT1B$Shp1 <- ifelse( DT1B$adj.P.Val < alpha & DT1B$logFC > folch, "T1", ifelse( DT1B$adj.P.Val < alpha & DT1B$logFC < -folch, "Bee","NS" ) )
	
	
	
#computing DE ASVs using fitZig
	physeq.i=subset_samples( physeq.I, Cond %in% c( "T2","Bee" ) )	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cond
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cond, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CondT2	-	CondBee, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable( fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames( DT.zig )
	DT.zig<-data.table( DT.zig, key="rn" )

	DT2B=DT.taxa[ DT.zig, ]

	DT2B$Col <- ifelse ( DT2B$adj.P.Val < alpha & abs(DT2B$logFC) > folch, DT2B$ColPhylum,"NS" )
	LIST.2B=DT2B[ DT2B$Col != "NS" ]$rn
#	DT2B$Shp1 <- ifelse( DT2B$adj.P.Val < alpha & DT2B$logFC > folch, "T2", ifelse( DT2B$adj.P.Val < alpha & DT2B$logFC < -folch, "Bee","NS" ) )
	


# transformt counts to rel. abundance
	physeq.ra = transform_sample_counts( physeq.I, function( x ) x/sum( x ) )

# computing arithmetic mean of rel. abund. by genotype
	DT.ra=data.table( otu_table( physeq.ra ), keep.rownames=T ) 
	DT.melt=data.table( melt( otu_table( physeq.ra ) ), keep.rownames=F )
	DT_sd=data.table( data.frame( sample_data( physeq.ra ) ), keep.rownames=T, key="rn" )
	DT.merge=merge( DT.melt, DT_sd, by.x="Var2", by.y="rn" )

	list.T1<-c( unique( DT.merge[ DT.merge$Cond=="T1", ]$Var2 ) )
	DT.ra$T1.mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.T1 ]
	
	list.T2<-c( unique( DT.merge[ DT.merge$Cond=="T2", ]$Var2 ) )
	DT.ra$T2.mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.T2 ]
	
	list.Bee<-c( unique( DT.merge[ DT.merge$Cond=="Bee", ]$Var2 ) )
	DT.ra$Bee.mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.Bee ]
	
	list.all<-c( unique( DT.merge[ DT.merge$Cond %in% c( "T1","T2","Bee" ) ]$Var2 ) )
	DT.ra$mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.all ]
	
	DT.mean=DT.ra[,c("rn","T1.mean","T2.mean","Bee.mean","mean")]
	DT.mean=merge(DT.mean, DT.taxa, by.x="rn", by.y="rn") 
	DT.mean=DT.mean[DT.mean$mean!=0]
#	DT.mean=DT.mean[DT.mean$mean>0.00001]

	LIST.SIG<-unique( c( LIST.12,LIST.1B,LIST.2B ) )
	

	DT.mean$lab=ifelse( DT.mean$rn %in% LIST.SIG, DT.mean$Genus,"" )
	DT.mean$lab=ifelse( DT.mean$mean > 0.01, DT.mean$Genus,"" )
	DT.mean$Col=ifelse( DT.mean$rn %in% LIST.SIG, DT.mean$ColPhylum,"NS" )

	minZ=1
	maxZ=10
	RAN=c( minZ,maxZ )
	minX=1e-07
	maxX=0.80
	LIM=c( minX,maxX )
	SEQ=seq( 0.05,0.8,0.05 )
	
#	DT.mean=DT.mean[DT.mean$mean > 0.05]

	pp1=ggtern( DT.mean,mapping=aes( T1.mean,T2.mean,Bee.mean,colour=Col ) )
	p1=pp1 + geom_point( aes(size=mean),alpha=.8 )+ theme_bw ( ) + 
	scale_color_manual( values=Palette.phylum )+ 
	scale_size_continuous( range=RAN,breaks=SEQ ) + 
	geom_text( aes( label=lab ) )+
	theme_change + add_seg + labs( x="T1",y="T2",z="Bee")

#	pdf("ternary_plot_all_days.pdf", useDingbats=FALSE)
	print(p1)
#	dev.off()

	


























