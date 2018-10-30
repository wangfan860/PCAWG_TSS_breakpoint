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

#define pdf file name
file_name <- paste("/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/", deparse(substitute(args[1])),deparse(substitute(args[5])), ".pdf", sep="")
pdf(file_name, width = 16, height = 8)

#cnv data
data <- read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/cnv_8990_for_cnvheatmap.csv', header=T)
data1 <- data[data$Chromosome==chrom,]
data2 <- data1[data1$project_id==cancer_type,]
data3 <- subset(data2, select = -c(X, file_name,project_id, Segment_Mean))
colnames(data3) <- c("chromosome", "start", "end", "sample","segmean")

#expression data
expr <- read.csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/102418_normalized_tcga_expr_forheatmap.csv', header=T)
expr1 <- expr[expr$project_id==cancer_type,]
expr2 <- subset(expr1, select=c(barcode, project_id,gene))

#merge cnv and expr
merge <- merge(expr2, data3, by.x = "barcode", by.y = "barcode")
merge1 <- merge[order(-gene),]
merge1$sample <- seq.int(nrow(merge1))

#edit theme for cnvheatmap
new_theme <- theme(panel.spacing.y = unit(0, "lines"), panel.border = element_rect(colour = "white", fill=NA, size=0), panel.grid.minor=element_line(colour="white", size=0),panel.grid.major =element_line(colour="white", size=0))
genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg38",], max)
genomeBoundaries$chromStart <- 0
colnames(genomeBoundaries) <- c("chromosome", "end", "start")
genomeBoundaries$start <- start
genomeBoundaries$end <-  end
cnvheatmap <- cnSpec(merge1, y=genomeBoundaries[genomeBoundaries$chromosome==chrom,], plotLayer = list(new_theme), x_title_size=0, facet_lab_size = 5, CNscale='absolute' )

#Define layout for the plots (2 rows, 2 columns)
layt<-grid.layout(nrow=2,ncol=2,heights=c(8/10,2/10),widths=c(1/6,5/6),default.units=c('null','null'))
#Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout=layt))
print(py,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(cnvheatmap,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(px,vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()
