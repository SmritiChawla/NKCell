##loading libraries
library(Seurat)
library(viridis)
library(RColorBrewer)
library(pheatmap)

#####Processing single cell gene expression data from run1 
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

#####Processing single cell gene expression data from run2
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
colnames(mt1) = c("chip","Run","Selection","Tumor","NK","Final")
colnames(mt2) = c("chip","Run","Selection","Tumor","NK","Final")

m1 = as.matrix(paste(mt1$Selection,mt1$Final,sep="_"))
rownames(m1) = rownames(m1)
pos1 = which(m1[,1]=="TU-NK_TU")
expression_matrix = expression_matrix[,-pos1]

m2 = as.matrix(paste(mt2$Selection,mt2$Final,sep="_"))
rownames(m2) = rownames(m2)
pos2 = which(m2[,1]=="TU-NK_TU")
expression_matrix1 = expression_matrix1[,-pos2]

###Running Seurat pipeline

#Breast_cancer_run1
df1 <- CreateSeuratObject(expression_matrix, project = "Breast_cancer_run_1", min.cells = 5)
df1@meta.data$stim <- "Run1"

#Breast_cancer_run2
df2 <- CreateSeuratObject(expression_matrix1, project = "Breast_cancer_run_2", min.cells = 5)
df2@meta.data$stim <- "Run2"


objects = list()
objects[[1]] = df1
objects[[2]] = df2
for (i in 1:length(objects)) {
  objects[[i]] <- NormalizeData(objects[[i]], verbose = FALSE)
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", 
                                       verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = objects, dims = 1:30,k.filter = 100)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(combined) <- "integrated"

###Visualization and Clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
exp = as.matrix(combined@assays$integrated@data)

###Computing correlation between gene expression profiles of doublets and distance of doublets
dis = read.table("Run1_Run2_distance_info.csv",sep=",",header=T,stringsAsFactors = F,row.names = 1)
rownames(dis)[72:113] = gsub('(.*)_\\w+', '\\1',rownames(dis)[72:113])
dis = dis[,3:ncol(dis)]
pos = which(colnames(exp) %in% rownames(dis))
exp = exp[,pos]
dis <- dis[colnames(exp),,drop=FALSE]
exp = as.matrix(t(exp))
dis = as.matrix((dis))
cor = cor(exp,dis,method="pearson")

###Selecting genes based on defined threshold
threshold <- 0.25
idx <- apply(abs(cor) > threshold, 1, any)
correlate= cor[idx, ]

##Apply cutree to get several clusters                     
p=pheatmap(correlate,fontsize_row = 4.3,cluster_rows = T,fontsize_col = 8,show_colnames = T,cluster_cols = F,angle_col = 45)
cv=data.frame((cutree(p$tree_row, k=4)))
colnames(cv) = "Module"
cv = as.matrix(cv)
cv[cv[,1]==1] = "Module1"
cv[cv[,1]==2] = "Module2"
cv[cv[,1]==3] = "Module3"
cv[cv[,1]==4] = "Module4"
cv = data.frame(cv)

n=brewer.pal(n = 4, name = "Set1")
ann_col = list(
  Module = c(Module1 = "#377EB8", Module2 = "#4DAF4A",Module3="#E41A1C",Module4="#984EA3")
)
q_colors =  15
v_colors =  viridis(q_colors)
                    
##Plotting heatmap
pheatmap(correlate,fontsize_row = 4.2,cluster_rows = T,fontsize_col = 6,cellheight = 4,cellwidth = 10,show_colnames = T,cluster_cols = F,angle_col = 45,annotation_row = cv,annotation_colors = ann_col,color=colorRampPalette(v_colors)(100),fontsize = 8)


