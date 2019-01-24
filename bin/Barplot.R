#load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(ecodist)
library(rgl)

#load data
smooth_class <- read.delim("Smooth_class_lvl.txt", check.names=FALSE, header=T, stringsAsFactors=FALSE)
pustular_class <- read.delim("Pustular_class_lvl.txt", check.names=FALSE, header=T, stringsAsFactors=FALSE)

#melt data
mm <- melt(smooth_class)
mm1 <- melt(pustular_class)

#Reorder y-axis surface to deep 
mm$variable = with(mm, factor(variable, levels = rev(levels(variable)))) #smooth
mm1$variable = with(mm1, factor(variable, levels = rev(levels(variable)))) #pustular

#ggplot2 stacked barplot (Smooth)
t <- ggplot(mm, aes(x=variable, y=value, fill=factor(Class))) 
t <- t + geom_bar(stat="identity")
t <- t + coord_flip() + theme_bw()
t <- t + labs(title="", x="Mat Layer Depth (mm)", y="Relative Abundance (%)")
t <- t + ggtitle("A") + theme(plot.title = element_text(hjust = 0, size=50))
t <- t + scale_fill_hue(l=55, c=105)
t <- t + guides(fill=guide_legend(title="Class"))
t #plot

#ggplot2 stacked barplot (Pustular)
t1 <- ggplot(mm1, aes(x=variable, y=value, fill=factor(Class))) 
t1 <- t1 + geom_bar(stat="identity")
t1 <- t1 + coord_flip() + theme_bw()
t1 <- t1 + labs(title="", x="Mat Layer Depth (mm)", y="Relative Abundance (%)")
t1 <- t1 + ggtitle("B") + theme(plot.title = element_text(hjust = 0, size=50))
t1 <- t1 + scale_fill_hue(l=55, c=105)
t1 <- t1 + guides(fill=guide_legend(title="Class"))
t1 #plot



