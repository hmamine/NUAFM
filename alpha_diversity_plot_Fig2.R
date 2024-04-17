#***require packages for the analysis the analysis
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")

color_palette<-c("#00b2b2","#b02b76","#ffa500","#155832","#4f7b95")

set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.ticks.x=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
#	axis.text.x = element_text(angle=90, vjust=1, size=10, color="black"),
	axis.text.x=element_blank(),
#	axis.text.y=element_blank(),
	)
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
}
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
#Rarefication to even depth based on smallest sample size
rf=min(sample_sums(physeq.subset))
physeq.rrf=rarefy_even_depth(physeq.subset, rf, replace=TRUE, rngseed = 131)

#Ploting species richness	
	p_nat=plot_richness(physeq.rrf,"Concat","Cond" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	Ord=c("T1D1","T2D1","T1D3","T2D3","T1D5","T2D5","Bee")
	levels(p_nat$data$Concat)<-Ord
	p_nat$data$Concat <- ordered(p_nat$data$Concat, levels= Ord)
	
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat, y=value))+
	geom_point(size=2)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

subp0=ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=2)+
	scale_fill_manual(values=color_palette) +  
	ylab("Rhichness")
	
subp2=ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=2)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette)
	legend <- g_legend(subp0) 
	
#	DT2=subp2$data
#	kruskal.test(data=T1, value ~ Concat)
#	kwAllPairsConoverTest(value~as.factor(Concat),data=DT1, p.adjust.method="BH")
#	aov.out<-aov(data=DT2, value ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	
#	T2=subp2$data[subp2$data$Cond =="T2",]
#	kruskal.test(data=T2, value ~ Day)
#	kwAllPairsConoverTest(value~as.factor(Concat),data=T2, p.adjust.method="BH")
	
#	T3=subp2$data[subp2$data$Cond =="T3",]
#	kruskal.test(data=T3, value ~ Day)
#	kwAllPairsConoverTest(value~as.factor(Concat),data=T3, p.adjust.method="BH")



subp1=ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=2)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) + 
	ylab("Shannon index")

#	DT1=subp1$data
#	kruskal.test(data=DT1, value ~ Concat)


otu.mat.rrf=as( otu_table( physeq.rrf), "matrix" )
Faith.PD<-data.table(pd(t(otu.mat.rrf), TREE), keep.rownames=T, key="rn")
setnames(Faith.PD, "rn", "samples")
DT_PD<-merge(subp2$data[,c(-7:-8)], Faith.PD)
DT_PD$variable <- "FaithPD"

subp3 = ggplot(data=DT_PD, aes(x=Concat, y=PD))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=2)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) + 
	ylab("Phylogenetic diversity")


#	DT3=subp3$data
#	kruskal.test(data=DT3, PD ~ Concat)
#	kwAllPairsConoverTest(PD~as.factor(Concat),data=DT3, p.adjust.method="BH")
#	aov.out<-aov(data=DT3, PD ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	
#	T2=subp2$data[subp2$data$Cond =="T2",]
#	kruskal.test(data=T2, value ~ Day)
#	kwAllPairsConoverTest(value~as.factor(Concat),data=T2, p.adjust.method="BH")
	
#	T3=subp2$data[subp2$data$Cond =="T3",]
#	kruskal.test(data=T3, value ~ Day)
#	kwAllPairsConoverTest(value~as.factor(Concat),data=T3, p.adjust.method="BH")

#	pdf("alpha_diversity_day.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(subp1, subp2,subp3,legend, nrow=2, ncol=2)
#	dev.off()
           
            
         

