##Loading libaries
library(Seurat)
library(cowplot)
library(pheatmap)


##Processing single cell gene expression data from run1 
data <- read.csv("Breast_cancer_run1.csv",sep=",",header = T,stringsAsFactors = F,row.names = 1)
expression_matrix = as.matrix(data[,5:ncol(data)])
gnames = as.matrix(data[,1])
rownames(expression_matrix) = gnames
expression_matrix = expression_matrix[which(rowSums(expression_matrix)>0),]
expression_matrix = expression_matrix[,-c(63,68,96)]
colnames(expression_matrix) = gsub('(.*)_\\w+', '\\1',colnames(expression_matrix))
meta1 = read.csv("Run1_cell_metadata.csv",sep=",",header=T,stringsAsFactors = F,row.names = 1,check.names = F,strip.white = T)
mt1 = meta1[colnames(expression_matrix),]
pos = which(mt1[,6]==0)
mt1 = mt1[-pos,]
expression_matrix = expression_matrix[,-pos]
genes = read.table("protein_coding_gene_list.txt",sep="\t")
pos1 = which(rownames(expression_matrix) %in% genes[,1])
expression_matrix = as.matrix(expression_matrix[pos1,])
exprs_Genes = apply(expression_matrix, 1, function(x) sum(x > 5)) >= 10
pos = which(colSums(expression_matrix)>2000)
expression_matrix = expression_matrix[exprs_Genes,pos]

##Processing single cell gene expression data from run2
data1 <- read.csv("Breast_cancer_run2.csv",sep=",",header = T,stringsAsFactors = F,row.names = 1,strip.white=T)
expression_matrix1 = as.matrix(data1[,5:ncol(data1)])
gnames1 = as.matrix(data1[,1])
rownames(expression_matrix1) = gnames1
expression_matrix1 = expression_matrix1[which(rowSums(expression_matrix1)>0),]
expression_matrix1 = expression_matrix1[,-c(32,49,121)]
colnames(expression_matrix1) = gsub('(.*)_\\w+', '\\1',colnames(expression_matrix1))
genes = read.table("protein_coding_gene_list.txt",sep="\t")
pos1 = which(rownames(expression_matrix1) %in% genes[,1])
expression_matrix1 = as.matrix(expression_matrix1[pos1,])
exprs_Genes = apply(expression_matrix1, 1, function(x) sum(x > 5)) >= 10
pos = which(colSums(expression_matrix1)>2000)
expression_matrix1 = expression_matrix1[exprs_Genes,pos]


meta2 = read.table("Run2_cell_metadata.csv",sep=",",header=T,stringsAsFactors = F,row.names = 1,check.names = F,strip.white = T)
rownames(meta2) = gsub('(.*)_\\w+', '\\1',rownames(meta2))

mt2 = meta2[colnames(expression_matrix1),]
pos = which(mt2[,6]==0)

mt2 = mt2[-pos,]
expression_matrix1 = expression_matrix1[,-pos]

colnames(mt1) =c("chip","Run","Selection","Tumor","NK","Final")
colnames(mt2) =c("chip","Run","Selection","Tumor","NK","Final")
m1 = as.matrix(paste(mt1$Selection,mt1$Final,sep="/"))
rownames(m1) = rownames(m1)
pos1 = which(m1[,1]=="Cancer/Cancer")
expression_matrix = expression_matrix[,pos1]

m2 = as.matrix(paste(mt2$Selection,mt2$Final,sep="/"))
rownames(m2) = rownames(m2)
pos2 = which(m2[,1]=="Cancer/Cancer")
expression_matrix1 = expression_matrix1[,pos2]

#Breast_cancer_run1
df1 <- CreateSeuratObject(expression_matrix, project = "Breast_cancer_run_1", min.cells = 5)
df1@meta.data$stim<- "Run1"

#Breast_cancer_run2
df2 <- CreateSeuratObject(expression_matrix1, project = "Breast_cancer_run_2", min.cells = 5)
df2@meta.data$stim <- "Run2"
set.seed(100)
objects = list()

##Integrating Run1 & Run2
objects[[1]] = df1
objects[[2]] = df2

for (i in 1:length(objects)) {
  objects[[i]] <- NormalizeData(objects[[i]], verbose = FALSE)
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", 
                                       verbose = FALSE)
}

immune.anchors <- FindIntegrationAnchors(object.list = objects, dims = 1:30,k.filter = 5)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30,k.weight = 5)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


# Visualization
set.seed(100)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = F,cols=c("red","blue"),pt.size=3,label.size = 10)
p2