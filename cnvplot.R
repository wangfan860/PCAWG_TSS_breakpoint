#! /usr/bin/env Rscript
args <- commandArgs(trailing = TRUE)
start <- args[1]
end <- args[2]
gene <- args[3]

install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("devtools",repos = "http://cran.us.r-project.org")
library(devtools)
install_github("griffithlab/GenVisR")
library('GenVisR')
require(ggplot2)
pdf("/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/test0.pdf", width = 16, height = 8)
data=read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/cnv_8990_for_cnvheatmap.csv', header=T)
data1=data[data$Chromosome==11,]
data2=data1[data1$project_id=='TCGA-ESCA']
new_theme <- theme(panel.spacing.y = unit(0, "lines"), panel.border = element_rect(colour = "white", fill=NA, size=0), panel.grid.minor=element_line(colour="white", size=0),panel.grid.major =element_line(colour="white", size=0), plot.margin = margin(0.85,1 ,4 , 3, "in"))
genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg38",], max)
genomeBoundaries$chromStart <- 0
colnames(genomeBoundaries) <- c("chromosome", "end", "start")
genomeBoundaries$start <- 2000000
genomeBoundaries$end <-   2350000
cnSpec(Data2, y=genomeBoundaries[genomeBoundaries$chromosome=="chr11",], plotLayer = list(new_theme, geom_vline(xintercept = c(2150342, 2155364,2155439,2170833), colour="seagreen4", size=0.5, linetype=2)), x_title_size=0, facet_lab_size = 5, CNscale='absolute' )
dev.off()
install.packages(“gplots”,repos = "http://cran.us.r-project.org")
library(gplots)

exampleData <- matrix(log2(rexp(1000)/rexp(1000)),nrow=200)
evar <- apply(exampleData,1,var)
mostVariable <- exampleData[evar>quantile(evar,0.75),]
pdf("/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/test.pdf", width = 16, height = 8)
heatmap.2(mostVariable,trace=”none”,col=greenred(10),key.title=start, key.xlab=end, key.ylab=gene)
dev.off()
