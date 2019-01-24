#Archaea Amplicons bubbleplot/Biplot SharkBay 
#Aug 11, 2015
#RAWIII

#load libraries
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(FactoMineR)
library(factoextra)

#load data
smooth_class <- read.delim("Smooth_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F)
pustular_class <- read.delim("Pustular_class_lvl.txt", check.names=FALSE, header=T, stringsAsFactors=F)
p_class <- read.delim("Pustular_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)
s_class <- read.delim("Smooth_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)

##################
####Biplot/PCA####
##################

#take sqrt and transponse data
p_classSQ  <- sqrt(p_class)
p_classT <- t(p_classSQ)
s_classSQ  <- sqrt(s_class)
s_classT <- t(s_classSQ)

#PCA of data
pca_s <- PCA(s_classT, graph = FALSE)
pca_p <- PCA(p_classT, graph = FALSE)

p <- fviz_pca_biplot(pca_p) + labs(title= "", x="PC1 (53.5%)", y="PC2 (23.4%)")
p <- p + ggtitle("B") + geom_text(size=10)
p <- p + theme_bw() + theme(axis.text=element_text(size=25),
                              legend.title=element_text(size=25),
                              strip.text.x=element_text(size=25),
                              legend.text=element_text(size=25),
                              axis.title=element_text(size=25),
                              plot.title = element_text(hjust = 0, size=50))
p <- p + geom_text(size=10)
p

p1 <- fviz_pca_biplot(pca_s) + labs(title= "", x="PC1 (37.9%)", y="PC2 (19.2%)")
p1 <- p1 + ggtitle("D") + theme(plot.title = element_text(hjust = 0, size=50))
p1 <- p1 + theme_bw() + theme(axis.text=element_text(size=25),
                            legend.title=element_text(size=25),
                            strip.text.x=element_text(size=25),
                            legend.text=element_text(size=25),
                            axis.title=element_text(size=25),
                            plot.title = element_text(hjust = 0, size=50))
p1

###################
####Bubble plot####
###################

#melt data to rearrange for ggplot2
mp <- melt(pustular_class)
ms <- melt(smooth_class)

#Reorder y-axis surface to deep for Stacked Barplots
ms$variable = with(ms, factor(variable, levels = rev(levels(variable)))) #smooth
mp$variable = with(mp, factor(variable, levels = rev(levels(variable)))) #pustula

p2 <- ggplot(mp, aes(x=variable, y=Class, color=variable)) + geom_point(aes(size=log10(value))) 
p2 <- p2 + theme_bw() + theme(axis.text=element_text(size=25),
                              legend.title=element_text(size=25),
                              strip.text.x=element_text(size=25),
                              legend.text=element_text(size=25),
                              axis.title=element_text(size=25),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p2 <- p2 + labs(title= "", x="Pustular layers", y="Archaeal Class")
p2 <- p2 + scale_colour_discrete(guide = FALSE)
p2 <- p2 + ggtitle("A") + theme(plot.title = element_text(hjust = 0, size=50))
p2 <- p2 + coord_flip() + guides(fill=guide_legend(title="Class"))
p2

p3 <- ggplot(ms, aes(x=variable, y=Class, color=variable)) + geom_point(aes(size=log10(value))) 
p3 <- p3 + theme_bw() + theme(axis.text=element_text(size=25),
                              legend.title=element_text(size=25),
                              strip.text.x=element_text(size=25),
                              legend.text=element_text(size=25),
                              axis.title=element_text(size=25),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p3 <- p3 + labs(title= "", x="Smooth layers", y="Archaeal Class")
p3 <- p3 + scale_colour_discrete(guide = FALSE)
p3 <- p3 + ggtitle("C") + theme(plot.title = element_text(hjust = 0, size=50))
p3 <- p3 + coord_flip() + guides(fill=guide_legend(title="Class"))
p3

grid.arrange(p2, p, p3, p1, ncol=2)
grid.arrange(t, p, t1, p1, ncol=2)