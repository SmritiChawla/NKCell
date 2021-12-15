##loading libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


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
m1 = as.matrix(paste(mt1$Selection,mt1$Final,sep="_"))
rownames(m1) = rownames(m1)
pos1 = which(m1[,1]=="TU-NK_TU")
mt1 = mt1[-pos1,]
expression_matrix = expression_matrix[,-pos1]
m2 = as.matrix(paste(mt2$Selection,mt2$Final,sep="_"))
rownames(m2) = rownames(m2)
pos2 = which(m2[,1]=="TU-NK_TU")
mt2 = mt2[-pos2,]
expression_matrix1 = expression_matrix1[,-pos2]


#Breast_cancer_run1
ctrl <- CreateSeuratObject(expression_matrix, project = "Breast_cancer_run_1", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"


#Breast_cancer_run2
stim <- CreateSeuratObject(expression_matrix1, project = "Breast_cancer_run_2", min.cells = 5)
stim@meta.data$stim <- "STIM"
set.seed(100)
objects = list()


##Integrating Run1 & Run2
objects[[1]] = ctrl
objects[[2]] = stim

for (i in 1:length(objects)) {
  objects[[i]] <- NormalizeData(objects[[i]], verbose = FALSE)
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", 
                                       verbose = FALSE)
}

immune.anchors <- FindIntegrationAnchors(object.list = objects, dims = 1:30,k.filter = 100)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)



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
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
mt1 = mt1[,-c(4,5)]
colnames(mt1) = c("Chip","Run", "Selection.step" ,"Final.decision")
mt2 = mt2[,-c(4,5)]
colnames(mt2) = c("Chip","Run", "Selection.step" ,"Final.decision")
labels = rbind(mt1,mt2)
um= p1$data[,1:2]
lb = labels[rownames(um),]
final_labels = paste(lb$Selection.step,lb$Final.decision,sep="/")
data = cbind.data.frame(um[,1],um[,2],final_labels)
colnames(data) = c("UMAP1","UMAP2","CellType")


n=brewer.pal(n = 4, name = "Set1")

g=ggscatter(data, x = "UMAP1", y = "UMAP2",
          color = "CellType", palette = n,legend="right",size=2.5) 
g+guides(color = guide_legend(override.aes = list(size=2.5)))




##PCA based visualization of cluster 1
exp <- cbind( expression_matrix[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ] ,
              expression_matrix1[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ])
meta = labels[colnames(exp),]


cl = read.table("cluster1.csv",sep=",",header=T,stringsAsFactors = F,row.names=1)
pos = which(colnames(exp) %in% cl[,1])

exp = log2(exp+1)
ex = exp[,pos]

lb = meta[colnames(ex),]
cell_metadata = as.matrix(paste(lb$Selection.step,lb$Final.decision,sep="_"))

###PCA

pc = prcomp(ex)

data = cbind.data.frame(pc$rotation[,1:2],cell_metadata)
colnames(data) = c("PC1","PC2","labels")



library(RColorBrewer)
n=brewer.pal(n = 4, name = "Set1")

library(ggpubr)
g=ggscatter(data, x = "PC1", y = "PC2",
            color = "labels", palette = n,legend="right",size=3) 
g+guides(color = guide_legend(override.aes = list(size=5))) + theme(axis.title=element_text(size=30),axis.text.x = element_text(size = 30),axis.text.y = element_text(size = 30))









