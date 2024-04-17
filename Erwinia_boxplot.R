#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "metagenomeSeq", "PMCMRplus","agricolae", "magrittr", "RColorBrewer", "seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#b02b76","#ffa500","#00b2b2","#155832","#4f7b95")
color_palette1<-c("#ffa500","#00b2b2")
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.y = element_blank(),
	axis.title.x = element_blank(),
	)
	
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
	physeq.subset <- subset_samples( physeq.prune, Cond != "T3" )
	
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

###transform counts to rel. abund.	
	physeq.RA = transform_sample_counts(physeq.subset, function(x) x/sum(x))

###agglomerate at the family level
	physeq.agg=tax_glom(physeq.RA, "Genus") 

#Agglomeration by sample type
#	physeq_agg=merge_samples(physeq.RA, "concat1") 
#	physeq.Prot=subset_taxa( physeq_agg, Phylum == "Proteobacteria" )


	m=as( otu_table( physeq.agg ), "matrix" )
	DT_m <- m %>% 
	melt(id.vars = "variable", value.name="relabund")
	DT.m=data.table(DT_m, keep.rownames=F, key="Var2") 
	
	DT.taxa=data.table( as(TAXA,"matrix"), keep.rownames=T, key="rn")
	DT.SD=data.table( as(SD,"data.frame"), keep.rownames=T, key="rn")
	
	DT=merge(DT.m, DT.taxa, by.x="Var1",by.y="rn")
	DT=merge(DT, DT.SD, by.x="Var2",by.y="rn")

	ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
	DT$ColPhylum <- ifelse(!DT$Phylum %in% ColPhylum, "Other", ifelse(DT$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT$Phylum == "Firmicutes", "Firmicutes", ifelse (DT$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))
	
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	Palette.Genus<-getPalette(40)

#	DT$Col<-ifelse(DT$relabund > 0.01, paste0(DT$Family,DT$Genus,DT$Var1), "Others")

	Ord1=c("Bee","T1D1","T2D1","T1D3","T2D3","T1D5","T2D5")

#	DT1=DT[DT$Genus == "Pseudomonas",]
#	pp1=ggplot(DT1, aes(x=Concat, y=relabund, color=Cond))
#	pp1$data$Concat <- ordered(pp1$data$Concat, levels=Ord1 )
#	p1=pp1+ geom_boxplot()+ theme_bw()+ theme_new +
#	scale_color_manual(values=color_palette) +
#	ylab("Relative Abundance of Pseudomonas")	

	DT1=DT[DT$Genus == "Erwinia",]
#	DT1=DT[DT$Family == "Erwiniaceae",]
	pp2=ggplot(DT1, aes(x=Concat, y=relabund, fill=Cond))
	pp2$data$Concat <- ordered(pp2$data$Concat, levels=Ord1 )
	p2=pp2+ geom_boxplot()+ geom_point(size=2, color="black")+
	theme_bw()+ theme_new + scale_fill_manual(values=color_palette) +
	ylab("Relative Abundance of Erwinia")	

	DT2=p2$data
	print(kruskal.test(data=DT2, relabund ~ Concat))
#	kwAllPairsConoverTest(value~as.factor(Concat),data=DT1, p.adjust.method="BH")
#	aov.out<-aov(data=DT2, relabund ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))


#	DT1=DT[DT$Col ==  "EnterobacteriaceaeND945184b6386c192c0066e0a98a154780",]
#	pp3=ggplot(DT1, aes(x=Concat, y=relabund, color=Cond))
#	pp3$data$Concat <- ordered(pp3$data$Concat, levels=Ord1 )
#	p3=pp3+ geom_boxplot()+ theme_bw()+ theme_new +
#	scale_color_manual(values=color_palette) +
#	ylab("Relative Abundance of Enterobacteriaceae")
	
	
	DT1=DT[DT$Genus ==  "Pantoea",]
	pp4=ggplot(DT1, aes(x=Concat, y=relabund, fill=Cond))
	pp4$data$Concat <- ordered(pp4$data$Concat, levels=Ord1 )
	p4=pp4+ geom_boxplot()+ geom_point(size=2, color="black")+
	theme_bw()+ theme_new + scale_fill_manual(values=color_palette) +
	ylab("Relative Abundance of Pantoea")
	
	DT1=DT[DT$Genus ==  "Pseudomonas",]
	pp5=ggplot(DT1, aes(x=Concat, y=relabund, fill=Cond))
	pp5$data$Concat <- ordered(pp5$data$Concat, levels=Ord1 )
	p5=pp5+ geom_boxplot()+ geom_point(size=2, color="black")+
	theme_bw()+ theme_new + scale_fill_manual(values=color_palette) +
	ylab("Relative Abundance of Pseudomonas")
	
	gridExtra::grid.arrange(p5,p4,p2,nrow=2,ncol=2)
break()

	tab1=read.table( "table_amsC.txt", sep="\t", row.names=1, header=T )
	DT1=data.table( tab1, keep.rownames=T, key="rn")
	DT1=DT1[DT1$Treat != "T3"]
	pp5=ggplot(DT1, aes(x=Concat, y=value, fill=Treat))
	p5=pp5+ geom_boxplot()+ geom_point(size=2, color="black")+
	theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette)
	
	DT5=p5$data
	print(kruskal.test(data=DT5, value ~ Concat))
	kwAllPairsConoverTest(value~as.factor(Concat),data=DT5, p.adjust.method="BH")
#	aov.out<-aov(data=DT5, value ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	
#	tab2=read.table( "table_relative_amsC.txt", sep="\t", row.names=1, header=T)
#	DT2=data.table( tab2, keep.rownames=T, key="rn")
#	DT2=DT2[DT2$Treat != "T3"]
#	pp6=ggplot(DT2, aes(x=Concat, y=-log(value), fill=Treat))
#	p6=pp6+ geom_boxplot()+ geom_point(size=2, color="black")+
#	theme_bw()+ theme_new +
#	scale_fill_manual(values=color_palette1)
	
#	pdf("Erwinia_plots.pdf", useDingbats=FALSE)
	gridExtra::grid.arrange(p2,p5,nrow=2,ncol=2)
#	dev.off()
	
	
	
	
	
	
	
	

