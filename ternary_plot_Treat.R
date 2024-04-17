pkg=c( "plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools","ggtern","PMCMRplus","venn", "metagenomeSeq","UpSetR", "venn","seqtime" )

lapply( pkg, library, character.only = TRUE )

rm( list = ls( )[ grep( "*.*", ls( ) ) ] )
set.seed( 131 )

Palette.phylum <-c( "#EF5656","#88888a","#2BB065" )
Palette.phylum <-c( "#EF5656","#47B3DA","#88888a","#2BB065" )
Palette.phylum <-c( "#88888a","#2BB065" )

ColPhylum<-c( "Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria" )
Shape.sign<-c( 21,16,15 )
lines <-data.frame( x=c( 0.5,0,0.5 ),y=c( 0.5,0.5,0 ),z=c( 0,0.5,0.5 ), xend=c( 1,1,1 )/3, yend=c( 1,1,1 )/3, zend=c( 1,1,1 )/3 )

minZ=1
maxZ=20
RAN=c( minZ,maxZ )

minX=1e-07
maxX=0.80
LIM=c( minX,maxX )

SEQ=seq( 0.05,0.8,0.15 )


add_seg<-geom_segment( data=lines, aes( x,y,z, xend=xend, yend=yend, zend=zend ),color="grey", size=1,linetype="dotdash" )

theme_change <- theme(
	legend.position="none",
	plot.background = element_blank( ),
	panel.grid.minor = element_blank( ),
	panel.grid.major = element_blank( ),
	)
DT<-vector( "list" )
DT.out<-list( )
RN<-list( )

#log2 fold change
	folch=2
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
	
	physeq.i <- prune_samples( sample_sums( physeq.MC )>=1000, physeq.MC )
	
	DT.taxa=data.table( as( tax_table( physeq.i ),"matrix" ), keep.rownames=T, key="rn" )
	DT.taxa$ColPhylum <- ifelse( !DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Actinobacteriota", "Actinobacteria", ifelse(DT.taxa$Phylum == "Bacteroidota","Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes","Firmicutes", "Proteobacteria" ) ) ) )
	
	
#trimming OTUs with low occurance
	otu.I=otu_table( physeq.i )
	taxa.I=tax_table( physeq.i )
	sd.I=sample_data( physeq.i )
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
	
	physeq.I=subset_samples( physeq.trim, Day == "D5" )
	
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

	DT13=DT.taxa[ DT.zig, ]

	DT13$Col <- ifelse( DT13$adj.P.Val < alpha, DT13$ColPhylum,"NS" )
	LIST.13=DT13[ DT13$Col != "NS" ]$rn
	
#computing DE ASVs using fitZig
	physeq.i=subset_samples( physeq.I, Cond %in% c( "T1","T3" ) )	
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
	c.mat1 = makeContrasts ( CondT1	-	CondT3, levels = finalMod1 )
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable( fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames( DT.zig )
	DT.zig<-data.table( DT.zig, key="rn" )

	DT15=DT.taxa[ DT.zig, ]

	DT15$Col <- ifelse ( DT15$adj.P.Val < alpha, DT15$ColPhylum,"NS" )
	LIST.15=DT15[ DT15$Col != "NS" ]$rn
	
#computing DE ASVs using fitZig
	physeq.i=subset_samples( physeq.I, Cond %in% c( "T2","T3" ) )	
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
	c.mat1 = makeContrasts ( CondT2	-	CondT3, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable( fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames( DT.zig )
	DT.zig<-data.table( DT.zig, key="rn" )

	DT35=DT.taxa[ DT.zig, ]

	DT35$Col <- ifelse ( DT35$adj.P.Val < alpha, DT35$ColPhylum,"NS" )
	LIST.35=DT35[ DT35$Col != "NS" ]$rn

	 
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
	
	list.T3<-c( unique( DT.merge[ DT.merge$Cond=="T3", ]$Var2 ) )
	DT.ra$T3.mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.T3 ]
	
	list.all<-c( unique( DT.merge[ DT.merge$Cond %in% c( "T1","T2","T3" ) ]$Var2 ) )
	DT.ra$mean<-DT.ra[ ,.( rowMeans( .SD,na.rm=TRUE ) ),.SDcols = list.all ]
	
	DT.mean=DT.ra[,c("rn","T1.mean","T2.mean","T3.mean","mean")]
	DT.mean=merge(DT.mean, DT.taxa, by.x="rn", by.y="rn") 
	DT.mean=DT.mean[DT.mean$mean!=0]
#	DT.mean=DT.mean[DT.mean$mean>0.00001]

	LIST.SIG<-unique( c( LIST.13,LIST.15,LIST.35 ) )

	DT.mean$lab=ifelse( DT.mean$rn %in% LIST.SIG, DT.mean$Genus,"" )
	DT.mean$Col=ifelse( DT.mean$rn %in% LIST.SIG, DT.mean$ColPhylum,"NS" )
	
	pp1=ggtern( DT.mean,mapping=aes( T1.mean,T2.mean,T3.mean,colour=Col ) )
	p1=pp1 + geom_point( size=2,alpha=.8 )+ theme_bw ( ) + 
	scale_color_manual( values=Palette.phylum )+ 
	scale_size_continuous( range=RAN,limits=LIM,breaks=SEQ ) + 
	geom_text( aes( label=lab ) )+
	theme_change + add_seg + labs( x="T1",y="T2",z="T3")

#	pdf("ternary_plot_D5.pdf", useDingbats=FALSE)
	print(p1)
#	dev.off()

	


























