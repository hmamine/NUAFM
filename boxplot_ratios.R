#require packages for the analysis the analysis
pkg=c("ggplot2", "data.table", "reshape2","PMCMRplus","agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#b02b76","#ffa500","#00b2b2","#155832","#4f7b95")

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank()
	)

	tab=read.table("ratio_table.txt", sep="\t", row.names=1, header=1)
	DT=data.table(tab, keep.rownames=T, key="rn")


	Ord1<-c("T1D1","T2D1","T1D3","T2D3","T1D5","T2D5")

	pp1<-ggplot(data=DT, aes(x=concat, y=log(value.16S), fill=treat, color=treat))

	levels(pp1$data$concat)<-Ord1
	pp1$data$concat <- ordered(pp1$data$concat, levels= Ord1)

	p1=pp1+stat_summary(fun="mean", geom="crossbar", size=.5, color="black")+
	geom_dotplot(binaxis = "y", stackdir = "center",binwidth = 0.1,show.legend=F)+
	scale_fill_manual(values=color_palette) + 
	scale_color_manual(values=color_palette) + 
	theme_bw() + theme_new
	
	
	pp2<-ggplot(data=DT, aes(x=concat, y=log(value.ITS), fill=treat, color=treat))

	levels(pp2$data$concat)<-Ord1
	pp2$data$concat <- ordered(pp2$data$concat, levels= Ord1)

	p2=pp2+stat_summary(fun="mean", geom="crossbar", size=.5, color="black")+
	geom_dotplot(binaxis = "y", stackdir = "center",binwidth = 0.07,show.legend=T)+
	scale_fill_manual(values=color_palette) + 
	scale_color_manual(values=color_palette) + 
	theme_bw() + theme_new



#	print(kruskal.test(data=DT, log(value.16S) ~ concat))
#	kwAllPairsConoverTest(log(value.16S) ~ as.factor(concat) , data=DT,  p.adjust.method="BH" )
#	aov.out<-aov(data=DT, log(value.16S) ~ concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="concat"))

#	print(kruskal.test(data=DT, log(value.ITS) ~ concat))
#	kwAllPairsConoverTest(log(value.ITS) ~ as.factor(concat) , data=DT,  p.adjust.method="BH" )
#	aov.out<-aov(data=DT, log(value.ITS) ~ concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="concat"))

	pp3<-ggplot(data=DT, aes(x=concat, y=log2(value.16S/value.ITS), fill=treat, color=treat))
	levels(pp3$data$concat)<-Ord1
	pp3$data$concat <- ordered(pp3$data$concat, levels= Ord1)

	p3=pp3+stat_summary(fun="mean", geom="crossbar", size=.5, color="black")+
	geom_dotplot(binaxis = "y", stackdir = "center",binwidth = 0.1,show.legend=T)+
	scale_fill_manual(values=color_palette) + 
	scale_color_manual(values=color_palette) + 
	theme_bw() + theme_new

#	pdf("dotplot_16S_ITS_2.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1,p2,p3,nrow=2, ncol=2)
#	dev.off()



#	print(kruskal.test(data=DT, log(value.16S/value.ITS) ~ concat))
#	kwAllPairsConoverTest(log(value.16S/value.ITS) ~ as.factor(concat) , data=DT,  p.adjust.method="BH" )
#	aov.out<-aov(data=DT, log(value1.6S/value.ITS) ~ concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="concat"))
















