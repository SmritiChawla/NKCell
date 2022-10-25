##Load libraries
library(pheatmap)

##Source function
source("Limma.R")

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

##Process metadata
colnames(mt1) = c("chip","Run","Selection","Tumor","NK","Final")
colnames(mt2) = c("chip","Run","Selection","Tumor","NK","Final")
labels = rbind((mt1),(mt2))

##Merge two runs
exp <- cbind( expression_matrix[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ] ,
              expression_matrix1[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ])

meta = labels[colnames(exp),]


##Load Cluster 1 and Cluster 2
Cluster1 = read.csv("cluster1.csv",row.names = 1)
Cluster2 = read.csv("cluster2.csv",row.names = 1)


##Subset Cluster 1 cells
##Subset Cluster 2 cells
data1 = exp[,which(colnames(exp) %in% rownames(Cluster1))]
data2 = exp[,which(colnames(exp) %in% rownames(Cluster2))]
count = cbind(data1,data2)
final_data = data.frame(count)
group = as.factor(c(rep("Cluster1",ncol(data1)),rep("Cluster2",ncol(data2))))

##Run limma voom
L = list(final_data,group)
names(L) <- c("count","condt")
voom = run_voomlimma(L)
top.table = voom$tt
mat = log2(count+1)
upgenes = top.table[top.table$logFC>0,]
downgenes = top.table[top.table$logFC<0,]

##Find top 20 genes  
upgenes_1 <- rownames(head(upgenes[ order( upgenes$logFC,decreasing = T), ], n=20))
downgenes_1 <- rownames(head(downgenes[ order( downgenes$logFC), ], n=20))
exp1 = mat[upgenes_1,]
exp2 = mat[downgenes_1,]
ex= rbind(exp1,exp2)

meta = rbind(mt1,mt2)
cells = meta[colnames(ex),]
labels = data.frame(colnames(ex),cells,group)
labels = labels[,-c(5,6)]
rownames(labels) = labels[,1]
labels = labels[,-1]
colnames(labels) = c("Chip","Run","TimePoint1","TimePoint16","Cluster")
labels$Run <- as.character(labels$Run)
labels$Run[labels$Run == "Run 1"] <- "Run1"
labels$Run[labels$Run == "Run 2"] <- "Run2"
labels$Run <- as.factor(labels$Run)
Chip = data.frame(labels[,1])
rownames(Chip) = rownames(labels)
colnames(Chip) = "Chip"
Run = data.frame(labels[,2])
colnames(Run) = "Run"
rownames(Run) = rownames(labels)
TimePoint1 = data.frame(labels[,3])
rownames(TimePoint1) = rownames(labels)
colnames(TimePoint1) = "TimePoint1"
TimePoint1$TimePoint1<-gsub("\\-","_",TimePoint1$TimePoint1)
TimePoint16 =  data.frame(labels[,4])
rownames(TimePoint16) = rownames(labels)
colnames(TimePoint16) = "TimePoint16"
TimePoint16$TimePoint16<-gsub("\\-","_",TimePoint16$TimePoint16)
Cluster = data.frame(labels[,5])
rownames(Cluster) = rownames(labels)
colnames(Cluster) = "Cluster"
labels$TimePoint1<-gsub("\\-","_",labels$TimePoint1)
labels$TimePoint16<-gsub("\\-","_",labels$TimePoint16)
ann_colors = list(
  Chip = c(chip1="red",chip2="green",chip3="blue",chip4="pink",chip5="orange",chip6="purple",chip8="yellow"),
  Run = c(Run1 = "red",Run2 = "blue"),
  TimePoint1 = c(Cancer="blue"),
  TimePoint16 = c(Cancer="blue"),
  Cluster = c(Cluster1="green",Cluster2="yellow")
)

##Plot heatmap
pheatmap(ex,annotation_col = labels,fontsize=7.5,cluster_rows = F,cluster_cols=F,show_colnames = F,annotation_colors = ann_colors,color=colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))(50))



