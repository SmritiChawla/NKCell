##loading libraries
library(readxl)
library(scatterD3)

excel_sheets("Run1_Run2_intensity_file.xlsx")
master1<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=1)
master2<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=2)
master3<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=3)
master4<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=4)
master5<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=5)
master6<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=6)
master7<-read_excel("Run1_Run2_intensity_file.xlsx", sheet=7)
View(head(master1))
master1[,c(2,3)]<-scale(master1[,c(2,3)])
master2[,c(2,3)]<-scale(master2[,c(2,3)])
master3[,c(2,3)]<-scale(master3[,c(2,3)])
master4[,c(2,3)]<-scale(master4[,c(2,3)])
master5[,c(2,3)]<-scale(master5[,c(2,3)])
master6[,c(2,3)]<-scale(master6[,c(2,3)])
master7[,c(2,3)]<-scale(master7[,c(2,3)])

final_list <- rbind(master1,master2,master3,master4,master5,master6,master7)

colnames(final_list) = c("cells","NK","Tumor","Annotation","Final.Annotation")

labels = paste(final_list$Annotation,final_list$Final.Annotation,sep="/")
pos1 = which(labels=="0/0")
pos2= which(labels=="NK/0")
pos3 = which(labels=="NK/TU-NK")
pos4 = which(labels=="TU-NK/TU")
pos5 = which(labels=="TU/0")
final_list = final_list[-c(pos1,pos2,pos3,pos4,pos5),]


final_list=data.frame(final_list)
scatterD3(data = final_list, x = NK, y = Tumor,labels_size = 0, point_size = 22,
          col_var = Final.Annotation,lab=cells,
          xlab = "Intensity of NK channel", ylab = "Intensity of tumor channel",
          col_lab = "Final annotation",
          symbol_lab = "Manual transmission")

