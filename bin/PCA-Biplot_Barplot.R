#Archaea Amplicons bubbleplot/Biplot SharkBay 
#Aug 12, 2015
#RAWIII

#load libraries
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(FactoMineR)
library(factoextra)

#load data
smooth_class <- read.delim("Smooth_class_lvl_other2.txt", check.names=F, header=T, stringsAsFactors=F)
pustular_class <- read.delim("Pustular_class_lvl_other.txt", check.names=FALSE, header=T, stringsAsFactors=F)
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

#melt data
mm <- melt(smooth_class)
mm1 <- melt(pustular_class)

#Reorder y-axis surface to deep for Stacked Barplots
mm$variable = with(mm, factor(variable, levels = rev(levels(variable)))) #smooth
mm1$variable = with(mm1, factor(variable, levels = rev(levels(variable)))) #pustular

#ggplot2 stacked barplot (Smooth)
t <- ggplot(mm, aes(x=variable, y=value, fill=factor(Class))) 
t <- t + geom_bar(stat="identity")
t <- t + coord_flip() + theme_bw() +
  theme(axis.text=element_text(size=25),
        legend.title=element_text(size=25),
        strip.text.x=element_text(size=25),
        legend.text=element_text(size=25),
        axis.title=element_text(size=25))
t <- t + labs(title="", x="Mat Layer Depth (mm)", y="Relative Abundance (%)")
t <- t + ggtitle("A") + theme(plot.title = element_text(hjust = 0, size=50))
t <- t + scale_fill_hue(l=55, c=105)
t <- t + guides(fill=guide_legend(title="lover lady"))
t #plot

#ggplot2 stacked barplot (Pustular)
t1 <- ggplot(mm1, aes(x=variable, y=value, fill=factor(Class))) 
t1 <- t1 + geom_bar(stat="identity")
t1 <- t1 + coord_flip() + theme_bw() +
  theme(axis.text=element_text(size=25),
        legend.title=element_text(size=25),
        strip.text.x=element_text(size=25),
        legend.text=element_text(size=25),
        axis.title=element_text(size=25))
t1 <- t1 + labs(title="", x="Mat Layer Depth (mm)", y="Relative Abundance (%)")
t1 <- t1 + ggtitle("C") + theme(plot.title = element_text(hjust = 0, size=50))
t1 <- t1 + scale_fill_hue(l=55, c=105)
t1 <- t1 + guides(fill=guide_legend(title="Class"))
t1 #plot

grid.arrange(t, t1, ncol=2)

png("test.png", width = 15, height = 15, units = 'in', res = 300)
grid.arrange(t, t1, ncol=2)
dev.off()

grid.arrange(t, p, t1, p1, ncol=2)