##loading libraries
library(tidyverse)
library(ggplot2)

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

colnames(mt1) = c("chip","Run","Selection","Tumor","NK","Final")
colnames(mt2) = c("chip","Run","Selection","Tumor","NK","Final")
m1 = as.matrix(paste(mt1$Selection,mt1$Final,sep="/"))
rownames(m1) = rownames(m1)
pos1 = which(m1[,1]=="Cancer-NK/Cancer")
mt1 = mt1[-pos1,]

expression_matrix = expression_matrix[,-pos1]
m2 = as.matrix(paste(mt2$Selection,mt2$Final,sep="/"))
rownames(m2) = rownames(m2)
pos2 = which(m2[,1]=="Cancer-NK/Cancer")
mt2 = mt2[-pos2,]
expression_matrix1 = expression_matrix1[,-pos2]
exp <- cbind( expression_matrix[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ] ,
              expression_matrix1[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ])
labels = rbind(mt1,mt2)
meta = labels[colnames(exp),]
cell_metadata = as.matrix(paste(meta$Selection,meta$Final,sep="/"))
rownames(cell_metadata) = rownames(meta)
colnames(cell_metadata) = "Cell_labels"
exp = data.frame(exp)
exp = log2(exp+1)
ex = cbind.data.frame((t(exp)),cell_metadata)
colnames(ex)[8908] = "cell_type" 

##Computing correlation for three Ligand/Protein pairs

##ANXA1/EGFR
ex1 = ex[,c("ANXA1","EGFR","cell_type")]

##Subsetting different cells
pos = which(ex1[,3]=="NK/NK")
nk = ex1[pos,-3]
cor1 = cor(nk)

pos = which(ex1[,3]=="Cancer-NK/Cancer-NK")
db = ex1[pos,-3]
cor2 = cor(db)

pos = which(ex1[,3]=="Cancer/Cancer")
Tu = ex1[pos,-3]
cor3 = cor(Tu)

pos = which(ex1[,3]=="Cancer-NK/NK")
tunk = ex1[pos,-3]
cor4 = cor(tunk)

c1 = data.frame(CellType=c("NK/NK", "Cancer-NK/Cancer-NK", "Cancer/Cancer","Cancer-NK/NK"), 
                ANXA1_EGFR=c(cor1[1,2],cor2[1,2],cor3[1,2], cor4[1,2]) )


##HSP90AA1/EGFR
ex1 = ex[,c("HSP90AA1","EGFR","cell_type")]

##Subsetting different cells
pos = which(ex1[,3]=="NK/NK")
nk = ex1[pos,-3]
cor1 = cor(nk)

pos = which(ex1[,3]=="Cancer-NK/Cancer-NK")
db = ex1[pos,-3]
cor2 = cor(db)

pos = which(ex1[,3]=="Cancer/Cancer")
Tu = ex1[pos,-3]
cor3 = cor(Tu)

pos = which(ex1[,3]=="Cancer-NK/NK")
tunk = ex1[pos,-3]
cor4 = cor(tunk)


c2 <- data.frame(CellType=c("NK/NK", "Cancer-NK/Cancer-NK", "Cancer/Cancer","Cancer-NK/NK"), 
                 HSP90AA1_EGFR=c(cor1[1,2],cor2[1,2],cor3[1,2], cor4[1,2]) )


###CD24/SIGLEC10
ex1 = ex[,c("CD24","SIGLEC10","cell_type")]

##Subsetting different cells
pos = which(ex1[,3]=="NK/NK")
nk = ex1[pos,-3]
cor1 = cor(nk)

pos = which(ex1[,3]=="Cancer-NK/Cancer-NK")
db = ex1[pos,-3]
cor2 = cor(db)

pos = which(ex1[,3]=="Cancer/Cancer")
Tu = ex1[pos,-3]
cor3 = cor(Tu)

pos = which(ex1[,3]=="Cancer-NK/NK")
tunk = ex1[pos,-3]
cor4 = cor(tunk)
c3 <- data.frame(CellType=c("NK/NK", "Cancer-NK/Cancer-NK", "Cancer/Cancer","Cancer-NK/NK"), 
                 CD24_SIGLEC10=c(cor1[1,2],cor2[1,2],cor3[1,2], cor4[1,2]) )

final = cbind.data.frame(c1,c2,c3)
final = final[,-c(3,5)]


##Plotting barplots               
df_long <- reshape2::melt(final)
df_long$CellType = factor(df_long$CellType,levels=c("NK/NK","Cancer-NK/NK","Cancer-NK/Cancer-NK", "Cancer/Cancer"))
ggplot(df_long, aes(variable, value, fill = CellType)) + 
  geom_bar(stat="identity", position = "dodge",width = 0.5) + 
  scale_fill_manual(values=c(c("NK/NK"="#E41A1C","Cancer-NK/NK"="#377EB8","Cancer-NK/Cancer-NK"="#FFD92F","Cancer/Cancer"="#984EA3")))+theme_classic(base_size = 20) + theme(axis.title=element_text(size=10),axis.text.x = element_text(size =10,hjust=1),axis.text.y = element_text(size = 10)) +theme(legend.position="right")  +  theme(legend.text=element_text(size=10))+ylab("Correlation")

