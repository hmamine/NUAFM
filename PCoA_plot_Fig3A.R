#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMRplus","vegan","agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#00b2b2","#b02b76","#ffa500","#155832","#4f7b95")
color.2<-c("#ffa500","#00b2b2","#155832","#4f7b95")
shape.1=c(13,17,15,16)
shape.2=c(17,15,16)

theme_new <- theme (
	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)
	
theme_new1 <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	axis.title.y=element_blank(),
#	axis.title.x=element_blank(),
	legend.position="none",
#	axis.text.x = element_text(angle=90, vjust=1, size=8),
#	axis.text.x=element_blank(),
#	axis.text.y=element_blank(),
	)
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")


set.seed(123)

# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	TREE=read.tree("tree.nwk")

	physeq=phyloseq(OTU, TAXA, SD) 
	
	physeq.MC=subset_taxa( physeq, !Genus %in% c("Mitochondria","Chloroplast") )
	
	physeq.prune <- prune_samples(sample_sums(physeq.MC)>=1000, physeq.MC)
	physeq.subset <- subset_samples(physeq.prune, Cond != "T3")
	Taxa=tax_table(physeq.subset)
	Sd=sample_data(physeq.subset)
	
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.subset ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), Taxa, Sd ,TREE )

#computing unconstrained PCoA using Bray-Curtis distances

	dist_BC=distance( physeq_norm, dist1 )

	pcoa_BC=ordinate( physeq_norm, meth1 , dist_BC )

	pBC<-plot_ordination( physeq_norm, pcoa_BC, color="Cond", shape="Day")
	
	pBC$layers<-pBC$layers[-1]
	
	p1=pBC+geom_point(size=5)+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)
	legend <- g_legend(p1) 
	
	p1=pBC+geom_point(size=4)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)
	
#DT.all=data.table( as( sample_data( physeq_norm ), "data.frame" ),keep.rownames=T, key="rn" )
#with( DT.all, adonis ( dist_BC ~ Comp ) )

#computing constrained PCoA	

	dist_BC2=distance( physeq_norm, dist1 )
	
	cpcoa_BC=ordinate( physeq_norm, meth2 , dist_BC2, ~Concat )

	cpBC<-plot_ordination( physeq_norm, cpcoa_BC, shape="Day", color="Cond")
	
	cpBC$layers<-cpBC$layers[-1]
	
	p2=cpBC+geom_point(size=4)+theme_bw()+theme_new1+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)
	
#	pdf("beta_plot.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p2,legend, nrow=2, ncol=2)
#	dev.off()

#	print(anova.cca(cpcoa_BC, by="axis", perm.max=1000))
	
	DT1=data.table(as(sample_data(physeq_norm),"data.frame"),keep.rownames=T, key="rn")
#	print(with(DT1, adonis2(dist_BC2 ~ Concat3)))

	physeq_norm2=subset_samples(physeq_norm, Cond !="Bee")
	dist_BC3=distance( physeq_norm2, dist1 )
	DT2=data.table(as(sample_data(physeq_norm2),"data.frame"),keep.rownames=T, key="rn")
	print(with(DT2, adonis2(dist_BC3 ~ Cond)))




