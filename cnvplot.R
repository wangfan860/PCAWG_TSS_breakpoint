#! /usr/bin/env Rscript
args <- commandArgs(trailing = TRUE)
cancer_type <- args[1]
chrom <- args[2]
start <- args[3]
end <- args[4]
gene <- args[5]

install.packages("grid",repos = "http://cran.us.r-project.org")
install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("devtools",repos = "http://cran.us.r-project.org")
install.packages("gplots",repos = "http://cran.us.r-project.org")
library(grid)
library(gplots)
library(devtools)
install_github("griffithlab/GenVisR")
library('GenVisR')
require(ggplot2)
source("http://bioconductor.org/biocLite.R")
biocLite("Sushi")
library('Sushi')
#define pdf file name
file_name <- paste("/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/", args[1],"_",args[5], ".pdf", sep="")
pdf(file_name, width = 8, height = 10)

#cnv data
data <- read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/cnv_8990_for_cnvheatmap.csv', header=T)
data1 <- data[data$Chromosome==chrom,]
data2 <- data1[data1$project_id==cancer_type,]
data3 <- subset(data2, select = -c(X, file_name,project_id, Segment_Mean))

#expression data
expr <- read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/tcga_expression_8990_top1000gene_for_cnvheatmap.csv', header=T)
expr1 <- expr[expr$project_id==cancer_type,]
expr2 <- subset(expr1, select=c('barcode', 'project_id',as.character(gene)))

#merge cnv and expr
merge <- merge(expr2, data3, by.x = "barcode", by.y = "barcode")
colnames(merge) <- c("barcode","project_id","gene","chromosome", "start", "end","segmean")
merge1 <- merge[order(-merge$gene),]
merge1$sample <- seq.int(nrow(merge1))

#cnvheatmap
new_theme <- theme(panel.spacing.y = unit(0, "lines"), panel.border = element_rect(colour = "white", fill=NA, size=0), panel.grid.minor=element_line(colour="white", size=0),panel.grid.major =element_line(colour="white", size=0))
genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg38",], max)
genomeBoundaries$chromStart <- 0
colnames(genomeBoundaries) <- c("chromosome", "end", "start")
genomeBoundaries$start <- start
genomeBoundaries$end <-  end
cnvheatmap <- cnSpec(merge1, y=genomeBoundaries[genomeBoundaries$chromosome==paste("chr",as.character(chrom),sep=""),], plotLayer = list(new_theme), x_title_size=0, facet_lab_size = 5, CNscale='absolute')

#expression_heatmap
data1=merge1[c('gene')]
data1$Gene=data1$gene
#data1_log=log(data1+1)
my_palette <- colorRampPalette(c("white",  "red"))(n = 601)
lmat = rbind(c(0,2), c(3,1), c(4,1),c(0,1))
lwid = c(1.5,1.5)
lhei = c(0.5,2.5,2.5,2.5)
expheatmap <- heatmap.2(as.matrix(unique(data1)),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,  Rowv=FALSE,    # use on color palette defined earlier
          Colv=FALSE,lmat = lmat, lwid = lwid, lhei = lhei)
#gene positions
gene=read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gencode_hg38_geneposition_forplot.txt', sep='\t', header=T)
chromstart = start
chromend   = end
pg = plotGenes(gene,chrom,chromstart,chromend,maxrows=1,bheight=0.15,
               plotgenetype="arrow",bentline=F,labeloffset=0.45 ,fontsize=0.6,
               arrowlength = 0,labeltext=T, addarrow=F,col='yellow')
labelgenome(chrom, chromstart,chromend,n=7,scale="Mb")
#Define layout for the plots (2 rows, 2 columns)
layt<-grid.layout(nrow=2,ncol=2,heights=c(8/10,2/10),widths=c(3/8,5/8),default.units=c('null','null'))
#Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout=layt))
print(expheatmap,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(cnvheatmap,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(pg,vp=viewport(layout.pos.row=2,layout.pos.col=2))

dev.off()
