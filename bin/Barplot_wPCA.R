#Archaea PCR/Barplots SharkBay 
#Aug 7, 2015
#RAWIII

#load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(ecodist)
library(rgl)
library(grid)
library(gridExtra)

#load data
smooth_class <- read.delim("Smooth_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F)
pustular_class <- read.delim("Pustular_class_lvl.txt", check.names=FALSE, header=T, stringsAsFactors=F)
genus <- read.delim("Genus_lvl_archaea.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)
s_genus <- genus[, -c(1,2,3,4,5)]
p_genus <- genus[, -c(6,7,8,9,10,11,12,13,14,15)]
p_class <- read.delim("Pustular_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)
p_class_rep <- read.delim("Pustular_mats_replicates.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)
s_class <- read.delim("Smooth_class_lvl.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)
s_class_rep <- read.delim("Smooth_mats_replicates.txt", check.names=F, header=T, stringsAsFactors=F, row.names=1)

  
  
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
t <- t + guides(fill=guide_legend(title="Class"))
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



library(FactoMineR)
library(factoextra)

p_classSQ  <- sqrt(p_class)
p_classT <- t(p_classSQ)

s_classSQ  <- sqrt(s_class)
s_classT <- t(s_classSQ)

p_classSQ_rep  <- sqrt(p_class_rep)
p_classTR <- t(p_classSQ_rep)

s_classSQ_rep  <- sqrt(s_class_rep)
s_classTR <- t(s_classSQ_rep)


res.pca_s <- PCA(s_classT, graph = FALSE)

res.pca <- PCA(p_class, graph = FALSE)
res.pca <- PCA(p_classT, graph = FALSE)
res.pca <- PCA(p_classT, graph = FALSE)
res.pca1 <- PCA(p_classTR, graph = FALSE)
res.pca2 <- PCA(s_classTR, graph = FALSE)
print(res.pca)


fviz_pca_biplot(res.pca_s) + theme_bw() + labs(title= "", x="PC1 (53.9%)", y="PC2 (24.9%)")
fviz_pca_biplot(res.pca2) + theme_bw() + labs(title= "", x="PC1 (53.9%)", y="PC2 (24.9%)")

fviz_pca_biplot(res.pca) + theme_bw() + labs(title= "", x="PC1 (53.9%)", y="PC2 (24.9%)")
fviz_pca_biplot(res.pca1) + theme_bw() + labs(title= "", x="PC1 (40.1%)", y="PC2 (16.3%)")

###############
###FactoMine###
###############

#variances
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

#scree plot
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")
fviz_screeplot(res.pca, ncp=10)

plot.PCA(x, axes = c(1,2), choix=c("ind", "var"))

plot(res.pca, choix = "var")
fviz_pca_var(res.pca)
fviz_pca_ind(res.pca)
y <- fviz_pca_biplot(res.pca),  geom = "")
fviz_pca_ind(res.pca, col.ind="cos2")

plot(res.pca, choix = "ind")
fviz_pca_biplot(res.pca,  geom = "text") + theme_bw()
fviz_pca_biplot(res.pca) + theme_bw()



fviz_pca_ind(res.pca,  col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.50)
fviz_pca_biplot(res.pca,  col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.50)

fviz_pca_ind(res.pca, habillage = 13,
             addEllipses =TRUE, ellipse.level = 0.68) +
  scale_color_brewer(palette="Dark2")
#########
###PCA###
#########

#Format data for PCA
data.hel  <- sqrt(p_genus)
data.hel.pca  <- rda(t(data.hel))
p <- length(data.hel.pca$CA$eig)
data.hel.pca.sc1 <- scores(data.hel.pca, display="wa", scaling=1, choices=c(1:p))
variance = (data.hel.pca$CA$eig / sum(data.hel.pca$CA$eig))*100

#cluster groups for PCA bray-curtis
source("pvclust_bcdist.R")
pathways_wide.hel.pv_fit <- pvclust(as.matrix(p_genus), method.hclust="ward", method.dist="bray–curtis", n=1)

# look at fit and decide cut height
plot(pathways_wide.hel.pv_fit) 
pathways_wide.hel.groups <- cutree(pathways_wide.hel.pv_fit$hclust, h=0.75) # slice dendrogram for groups

#plot PCA
p1 <- qplot(data.hel.pca.sc1[,1], 
            data.hel.pca.sc1[,2], 
            label=rownames(data.hel.pca.sc1), 
            size=2, geom=c("point"), 
            xlab= paste("PC1 (", round(variance[1],2) ," % Variance)"), 
            ylab= paste("PC2 (", round(variance[2],2) ," % Variance)"), 
            color= factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() + xlim(-15,15) + theme(legend.position="none") 
p1