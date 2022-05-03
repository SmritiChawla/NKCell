##load libraries
library(ggplot2)


##Load mRNA half life file
final = read.csv("K562_Amit_PR_half_life.csv",sep=",",header = T,stringsAsFactors = F,row.names = 1)


###Data processing
gM = read.csv("Module_genes.csv",sep=",")
pos = which(final$hgnc_symbol %in% gM[,1])
df1 = as.data.frame(final[pos,])
df1 = unique(df1)
df2 = as.data.frame(final[-pos,])
df2 = unique(df2)


colnames(df1)[2] = "Value"
colnames(df2)[2] = "Value"

group = c(rep("Module",nrow(df1)),rep("Non-Module",nrow(df2)))

mat = cbind.data.frame(rbind.data.frame(df1,df2),group)
mat$Value = log2(mat$Value)

#Plotting
ggplot(mat, aes(x=group, y=Value,col=group)) + 
  geom_boxplot() + theme_classic(base_size = 15) + scale_color_brewer(palette="Dark2")


##Statistical test for significance
t.test(log2(df1$Value),log2(df2$Value))



