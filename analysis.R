#-----------------------------------------------------------------------------------
# Compare RNA-seq data of dc-hiOL samples with public data from Jäkel et al. 2019 
# doi: 10.1038/s41586-019-0903-2, GEO accession number: GSE118257
#
# Christian Thomas
# c.thomas@uni-muenster.de                                                                 
# 
# 2020-12-08
#------------------------------------------------------------------------------------   

# Load libraries
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(preprocessCore)
library(DESeq2)
library(openxlsx)

# Load dc-hiOL data
hiOL_anno <- read.csv("./data/samples.txt",sep = ",")
hiOL<-read.table("./data/hiOL_salmon.counts")

# Load sn-RNAseq data from Jäkel et al. 2019 (GSE118257)
geoanno<-read.table("./data/GSE118257_anno.txt", header=TRUE, stringsAsFactors = FALSE)
geoexpr<-read.table("./data/GSE118257_expr.txt.gz", header=TRUE, row.names=1, stringsAsFactors = FALSE)

# Get intersection of gene names from both datasets
hiOL<-hiOL[intersect(rownames(hiOL),rownames(geoexpr)),]
geo<-geoexpr[intersect(rownames(hiOL),rownames(geoexpr)),]
colnames(geo)<-geoanno$Celltypes

# Calculate average expression values for each sn-RNAseq cluster (CellType)
geomean<-as.data.frame(t(geo))
geomean$Celltypes<-geoanno$Celltypes
geomean<-as.data.frame(geomean %>% group_by(Celltypes) %>% summarise_all(mean))

# Remove Endothelial cells, Macrophage and Vasc_smooth_muscle
geomean<-geomean[-c(4,5,7,22),]

# Tidy data
rownames(geomean)<-geomean[,1]
geomean<-geomean[,-1]

# Normalize sn-RNAseq count data
geomean<-t(geomean)
geomean<-round(geomean*100,0)
geomean_norm<-data.frame(normalize.quantiles(as.matrix(geomean)))
colnames(geomean_norm)<-colnames(geomean)
rownames(geomean_norm)<-rownames(geomean)

# Combine both datasets
comb<-cbind(hiOL, geomean_norm)

# Normalize combined data
comb_norm<-data.frame(normalize.quantiles(as.matrix(comb)))
comb_norm<-round(comb_norm,0)
colnames(comb_norm)<-colnames(comb)
colnames(comb_norm)[14]<-"Microglia"
rownames(comb_norm)<-rownames(comb)

# Import count data in DEseq2 and apply a variance stabilizing transformation (VST)
coldata<-data.frame(batch=c(rep("hiOL",9),rep("snRNA",18)))
rownames(coldata)<-colnames(comb_norm)
dds <- DESeqDataSetFromMatrix(comb_norm, coldata, ~batch)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

# Import marker genes
genesofinterest<-read.xlsx("./data/marker_genes.xlsx")

mat<-mat[match(genesofinterest$Gene,rownames(mat)),]
colnames(mat)<-c(paste0("dc-hiOL",seq(1,3,1)),
                  rep("Fib",3),rep("pOL",3),
                  colnames(mat[,10:length(mat[1,])]))

# Row annotation
row_annotation<-data.frame(Celltype=genesofinterest$Celltype)
rownames(row_annotation)<-rownames(mat)
celltype_cols<-brewer.pal(6,"Dark2")
names(celltype_cols)<-unique(genesofinterest$Celltype)
anno_cols<-list(Celltype=celltype_cols)

pheatmap(mat,
         color=colorRampPalette(c("blue","blue","white","yellow"))(255),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         annotation_row = row_annotation,
         annotation_colors = anno_cols,
         cluster_rows = FALSE
)

