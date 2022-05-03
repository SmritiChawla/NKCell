##loading libraries
library(limma)
library(edgeR)

##Processing single cell data 1
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

##Processing single cell data 2
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
labels = rbind((mt1),(mt2))

exp <- cbind( expression_matrix[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ] ,
              expression_matrix1[ intersect(rownames(expression_matrix), rownames(expression_matrix1)), ])

meta = labels[colnames(exp),]

cell_metadata = as.matrix(paste(meta$Selection,meta$Final,sep="_"))
rownames(cell_metadata) = rownames(meta)
colnames(cell_metadata) = "Cell_labels"

pos = which(cell_metadata[,1]=="TU-NK_TU")
cell_metadata = as.matrix(cell_metadata[-pos,])
exp = exp[,-pos]


##Running limma
nk_killing = read.table("NK_killing_16062020.txt",sep="\t",header = F)
touching = read.table("touching_16062020.txt",sep="\t",header = F)
nt = read.table("not_touching_16062020.txt",sep="\t",header = F)


data1 = exp[,which(colnames(exp) %in% nk_killing[,1])]
data2 = exp[,which(colnames(exp) %in% touching[,1])]
data3 = exp[,which(colnames(exp) %in% nt[,1])]


counts = cbind(data1,data2,data3)
group = c(rep("NK_Killing",10),rep("Not_killing",106))
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
mm <- model.matrix(~0 + group)
y <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
contrast.matrix <- makeContrasts(
  KN = groupNK_Killing - groupNot_killing,
  levels = colnames(coef(fit))
)
tmp <- contrasts.fit(fit, contrast.matrix)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, n = Inf)
