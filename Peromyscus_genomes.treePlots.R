   
#### Plot statistics of Peromyscus genomes on the Peromyscus tree
#### Based on multiple functions and examples from Liam Revell (http://blog.phytools.org/)

library(phytools)
library(RColorBrewer)

#### Read in data
data<-read.table("Peromyscus_assemblies_data.txt", sep="\t", header=T, row.names=1)
tree <- read.tree("Peromyscus_genomes.tre")
ape::nodelabels()
tree<-rotateNodes(tree, 8)

#### Create the palette
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(3, "Dark2")
brewer.pal(3, "Dark2")

#### Label tree by sequencing technology
fmode<-setNames(as.factor(data$sequencing),rownames(data))
plotTree(tree,ftype="i",offset=0.6,fsize=1.5)
FMODE<-to.matrix(fmode,levels(fmode))
FMODE<-FMODE[tree$tip.label,]

tiplabels(pie=FMODE,piecol=palette(c("#7570B3","#D95F02")),cex=0.4)
legend("bottomleft",levels(fmode),pch=21,pt.bg=palette(c("#7570B3","#D95F02")),pt.cex=2)
axis(3,at=seq(0,12,by=2), pos=4.5)

### Draw barplots of contig and scaffold N50 values

contig<-setNames(data$contig_N50,rownames(data))
obj_contig <- plotTree.barplot(tree,contig,
                       args.plotTree=list(fsize=0.7,ftype="i",offset=0.6,fsize=1.5),
                       args.barplot=list(width=0.2, xlim=c(0,12)))
                       
scaffold<-setNames(data$scaffold_N50,rownames(data))                  
obj_scaffold <- plotTree.barplot(tree,log(scaffold),
                       args.plotTree=list(fsize=0.7,ftype="i",offset=0.6,fsize=1.5),
                       args.barplot=list(width=0.2, xlim=c(0,5)))                   
                       