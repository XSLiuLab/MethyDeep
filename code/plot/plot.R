###########################  饼图
library(ggplot2)
b<-data.frame(c("Illumina Human Methylation 27","Illumina Human Methylation 450","WGBS CpG site signal"),
              c(2595,9736,47))
colnames(b)<-c("type","number")
myLabel = as.vector(b$type)   
myLabel = paste0(myLabel,"(",b$number,")")   
ggplot(b, aes(x = "", y = number, fill = type)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top")+ 
  scale_fill_discrete(breaks = b$type, labels = myLabel)+ theme(axis.text.x = element_blank(),panel.background = element_blank())
################### 条形图
library(ggsci)
library(ggprism)
b<-data.frame(c("A","B","C","D","E"),
              c(13884,444,153,34,3))
colnames(b)<-c("type","number")

ggplot(b,aes(x=type,y=number,fill=type))+geom_bar(stat = 'identity')+theme_prism()+scale_fill_npg()+
  geom_text(aes(label=number, y=number+0.05), position=position_dodge(0.9), vjust=0)+theme(legend.position = 'none')

################  条形图

b<-read.csv("F:/result/result_11_15/process_data/methylation_data.csv")
b<-b[,-1]
colnames(b)<-c("Cancertype","number")
color<-c(pal_nejm("default", alpha =0.7)(8),pal_npg("nrc", alpha =0.7)(10),pal_aaas("default", alpha =0.7)(10),pal_jama("default", alpha =0.7)(5))
p<-ggplot(b,aes(x=Cancertype,y=number,fill=Cancertype))+geom_bar(stat = 'identity')+theme_prism()+scale_fill_manual( values = color)+
  geom_text(aes(label=number, y=number+0.05), position=position_dodge(0.9), vjust=0)+theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=.9))
ggsave(p,file="F:/result/result_11_15/plot/fig1/甲基化数目.pdf",width = 30,height = 20)

b<-read.csv("F:/result/result_11_15/process_data/cancer_sample.csv")
b<-b[,-1]
colnames(b)<-c("Cancertype","number")

b<-b[order(b[,2],decreasing = T),]
b[,1]<-factor(b[,1],levels = b[,1])

color<-c(pal_nejm("default", alpha =0.7)(8),pal_npg("nrc", alpha =0.7)(10),pal_aaas("default", alpha =0.7)(10),pal_jama("default", alpha =0.7)(5))
ggplot(b,aes(x=Cancertype,y=number,fill=Cancertype))+geom_bar(stat = 'identity')+theme_prism()+scale_fill_manual( values = color)+
  geom_text(aes(label=number, y=number+0.05), position=position_dodge(0.9), vjust=0)+theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=.9))
ggsave(p,file="F:/result/result_11_15/plot/fig1/甲基化数目.pdf",width = 30,height = 20)

############################### 密度图
library(data.table)
library(stringr)

a<-list.files("F:/result/result_11_15/process_data/wilcox_differen_result/")

b<-str_split_fixed(a,"\\.",2)[,1]
c<-data.frame("1",1)
for(i in 1:length(a)){
  data<-fread(paste0("F:/result/result_11_15/process_data/wilcox_differen_result/",a[i]),data.table=F)
  #data<-fread(paste0("F:/result/result_11_15/process_data/origin_differen_result/",a[i]),data.table=F)
  #print(length(setdiff(data1[,1],data[,2])))
   c[i,1]<-b[i]
   c[i,2]<-nrow(data)
}
ggplot(c,aes(x=X1))+geom_density(fill="#67C2A3", color=c("#67C2A3"), alpha=0.8)+theme_classic()+xlab("number of DMS")

c<-c[order(c[,2],decreasing = T),]
ggplot(c[1:5,],aes(x=X.1.,y=X1,fill=X.1.))+geom_bar(stat = 'identity')+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10))+theme_prism()+theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))+
  xlab("")+ylab("number of DMS")+theme(legend.position = 'none')


all_data<-c()
b<-seq(1,435,30)
b<-c(b,435)
for(j in 1:15){
  all_data<-c()
  if(j!=15){
  m<-b[j]
  n<-b[j+1]-1
  }else{
    m<-b[j]
    n<-b[j+1]
  }
  for(i in m:n){
    data<-fread(paste0("F:/result/result_11_15/process_data/wilcox_differen_result/",a[i]),data.table=F)
    if(nrow(data)!=0){
      all_data<-c(data[,2],all_data)
    }
  }
  c<-as.data.frame(table(all_data))
  write.csv(c,file = paste0("F:/result/result_11_15/process_data/CPG_number/",j,".csv"),row.names = F)
}

library(dplyr)
a<-list.files("F:/result/result_11_15/process_data/CPG_number/")
all_data<-fread(paste0("F:/result/result_11_15/process_data/CPG_number/",a[1]),data.table = F)

for(i in 2:15){
  data<-fread(paste0("F:/result/result_11_15/process_data/CPG_number/",a[i]),data.table = F)
  all_data<-inner_join(all_data,data,by="all_data")
  
  }

a<-apply(all_data[,-1],1,sum)
a<-as.data.frame(a)
a[,2]<-all_data[,1]
ggplot(a,aes(x=a))+geom_density(fill="#FC8A61", color="#FC8A61", alpha=0.8)+theme_classic()+xlab("number of Recurrent DMS")
write.csv(a,file = "F:/result/result_11_15/supplement_data/wilcox_REcurren.csv",row.names = F)
###############################  champ
a<-list.files("F:/result/result_11_15/process_data/origin_differen_result/")

b<-str_split_fixed(a,"\\.",2)[,1]
c<-data.frame("1",1)
for(i in 1:length(a)){
  #data<-fread(paste0("F:/result/result_11_15/process_data/wilcox_differen_result/",a[i]),data.table=F)
  data<-fread(paste0("F:/result/result_11_15/process_data/origin_differen_result/",a[i]),data.table=F)
  #print(length(setdiff(data1[,1],data[,2])))
  c[i,1]<-b[i]
  c[i,2]<-nrow(data)
}
c1<-c[which(c[,2]<1000),]
ggplot(c1,aes(x=X1))+geom_density(fill="#67C2A3", color=c("#67C2A3"), alpha=0.8)+theme_classic()+xlab("number of DMS")

c<-c[order(c[,2],decreasing = T),]
ggplot(c[1:5,],aes(x=X.1.,y=X1,fill=X.1.))+geom_bar(stat = 'identity')+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10))+theme_prism()+theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))+
  xlab("")+ylab("number of DMS")+theme(legend.position = 'none')



all_data<-c()
for(i in 1:length(a)){
  #data<-fread(paste0("F:/result/result_11_15/process_data/wilcox_differen_result/",a[i]),data.table=F)
  data<-fread(paste0("F:/result/result_11_15/process_data/origin_differen_result/",a[i]),data.table=F)
  all_data<-c(all_data,data[,1])
}
c<-as.data.frame(table(all_data))
ggplot(c,aes(x=Freq))+geom_density(fill="#FC8A61", color=c("#FC8A61"), alpha=0.8)+theme_classic()+xlab("number of Recurrent DMS")
write.csv(c,file = "F:/result/result_11_15/supplement_data/champ_REcurren.csv",row.names = F)
############################### 热图
library(data.table)
library(stringr)

a<-list.files("F:/result/result_11_15/process_data/origin_differen_result/")
b<-str_split_fixed(a,"\\.",2)[,1]
b<-str_split_fixed(b,"_",2)[,1]
b<-unique(b)
b<-c(b,"TCGA-UVM")
data<-matrix(1,ncol=30,nrow=30)
colnames(data)<-rownames(data)<-b

library(corrplot)
library(ggprism)
corrplot(corr =data,method = "circle",type = "lower",col="red",diag=FALSE,tl.col="black")+theme_prism()


######################  保留100个

a<-list.files("F:/result/result_11_15/process_data/origin_differen_result/")
b<-str_split_fixed(a,"\\.",2)[,1]
c<-data.frame("1",1)

recurren_number<-read.csv("F:/result/result_11_15/supplement_data/wilcox_REcurren.csv")
recurren_number<-recurren_number[order(recurren_number[,1],decreasing = T),]

i=1
all_cg<-data.frame()
for(i in 1:length(a)){
  data<-fread(paste0("F:/result/result_11_15/process_data/origin_differen_result/",a[i]),data.table=F)
  if(nrow(data)>=100){
    data[,2]<-abs(data[,2])
    data<-data[,1:2]
    data<-data[order(data[,2],decreasing = T),]
    all_cg[i,1]<-b[i]
    all_cg[i,2:101]<-data[1:100,1]
  }else{
    wilcox_data<-fread(paste0("F:/result/result_11_15/process_data/wilcox_differen_result/",a[i]),data.table=F)
    colnames(wilcox_data)[2]<-"V2"
    wilcox_data<-left_join(recurren_number,wilcox_data,by="V2")
    wilcox_data<-na.omit(wilcox_data)
    cg<-c(data[,1],wilcox_data[1:100,2])
    cg<-cg[1:100]
    all_cg[i,1]<-b[i]
    all_cg[i,2:101]<-cg
  }
  i<-i+1
  
}
write.csv(all_cg,file = "F:/result/result_11_15/supplement_data/cg100.csv",row.names = F)

d<-c()
for(i in 2:101){
  d<-c(d,all_cg[,i])
}
d<-unique(d)
write.csv(d,file = "F:/result/result_11_15/supplement_data/all_cg.csv",row.names = F)

#################################  分组条形图
library(ggplot2)
library(ggprism)
library(ggsci)
a<-list.files("f:/result/result_11_15/process_data/ML_result_11_23/")
a<-a[grep("11_23.txt$",a)]
c<-str_split_fixed(a,"\\.",2)[,1]
d<-data.frame()
for(i in 1:length(a)){
  b<-read.table(paste0("f:/result/result_11_15/process_data/ML_result_11_23/",a[i]))
  b["method"]<-c[i]
  d<-rbind(d,b)
}
d[,1]<-str_split_fixed(d[,1],":",2)[,1]
d[which(d[,3]=="K-紧邻_11_23"),3]<-"KNN"
d[which(d[,3]=="决策树_11_23"),3]<-"Decision Tree"
d[which(d[,3]=="朴素贝叶斯_11_23"),3]<-"Naive Baye"
d[which(d[,3]=="随机森林_11_23"),3]<-"Random forest"
d[which(d[,3]=="XGBOOST_11_23"),3]<-"XGBOOST"
d[which(d[,3]=="ExtraTree_11_23"),3]<-"Extratree"

# d$method<-as.factor(d$method)
# 
# ggplot(d,aes(x=V1,y=V2,color=method,group=method))+
#   theme_prism()+geom_smooth(size=1.5)+geom_point(size=2)+ylab("")+xlab("Model evaluation index")

d[,3]<-factor(d[,3],levels = c("XGBOOST","Extratree","Random forest","KNN","Naive Baye","Decision Tree"))

ggplot(d,aes(x=method,y=V2,fill=V1))+geom_bar(stat = 'identity',position = position_dodge(0.9))+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10))+theme_prism()+theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))+
  xlab("")+ylab("number of DMS")+geom_hline(aes(yintercept=0.93),color="red")

+theme(legend.position = 'none') ## 去除legend
###############  figure4

a<-list.files("f:/result/result_11_15/process_data/ML_different_number_11_23/")
d<-data.frame()
for(i in c("ExtraTree","XGBOOST","随机森林")){
  for(j in c(30)){
    b<-read.table(paste0("f:/result/result_11_15/process_data/ML_different_number_11_23/",i,j,"_11_23.txt"))
    b["number"]<-j
    b["method"]<-i
    d<-rbind(d,b)
  }
}
d[d[,4]=="随机森林",4]<-"Random forest"

d1<-d
d1<-d1[,-3]

a<-list.files("f:/result/result_11_15/process_data/DNN_different_number_11_23/")
a<-a[grep("30_11_23.txt$",a)]
c<-str_split_fixed(a,"\\.",2)[,1]
d<-data.frame()
for(i in 1:length(a)){
  b<-read.table(paste0("f:/result/result_11_15/process_data/DNN_different_number_11_23/",a[i]))
  b["method"]<-c[i]
  d<-rbind(d,b)
}
d[,1]<-str_split_fixed(d[,1],":",2)[,1]
d[which(d[,3]=="DNN_Extra30_11_23"),3]<-"DNN_Extra"
d[which(d[,3]=="DNN_RF30_11_23"),3]<-"DNN_RF"
d[which(d[,3]=="DNN_Xgboost30_11_23"),3]<-"DNN_Xgboost"
d<-d[-c(5:8),]

d<-rbind(d,d1)

d[,3]<-factor(d[,3],levels = c("DNN_Xgboost","XGBOOST","DNN_Extra","ExtraTree", "DNN_RF", "Random forest"))

ggplot(d,aes(x=method,y=V2,fill=V1))+geom_bar(stat = 'identity',position = position_dodge(0.9))+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10))+theme_prism()+theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))+
  xlab("")+ylab("number of DMS")+geom_hline(aes(yintercept=0.93),color="red")

+theme(legend.position = 'none') ## 去除legend

############################## upset
library(UpSetR)
RF<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/RF_feature_11_23.csv")
EXtra<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Extratree_feature_11_23.csv")
XGBOOST<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Xgboost_feature_11_23.csv")


RF<-RF[1:100,]
EXtra<-EXtra[1:100,]
XGBOOST<-XGBOOST[1:100,]
upset_list<-list(RF[,1],EXtra[,1],XGBOOST[,1])
names(upset_list) <- c("Random forest","Extratree","XGBoost")
color<-c(pal_nejm("default", alpha =0.7)(8),pal_npg("nrc", alpha =0.7)(10),pal_aaas("default", alpha =0.7)(10),pal_jama("default", alpha =0.7)(5))

#作图
upset(fromList(upset_list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 ,# 文字标签的大小
      
      main.bar.color=color[1:6],
      sets.bar.color=color[7:9]
)

Reduce(intersect, list(RF[,1],EXtra[,1],XGBOOST[,1]))

############################## 折线图

a<-list.files("f:/result/result_11_15/process_data/ML_different_number_11_23/")
d<-data.frame()
for(i in c("ExtraTree","XGBOOST","随机森林")){
  for(j in c(5,10,15,20,25,30,35,40,45,50)){
    b<-read.table(paste0("f:/result/result_11_15/process_data/ML_different_number_11_23/",i,j,"_11_23.txt"))
    b["number"]<-j
    b["method"]<-i
    d<-rbind(d,b)
  }
}
d[d[,4]=="随机森林",4]<-"Random forest"
a<-unique(d[,1])
for(i in a){
  data<-d[d[,1]==i,]
  p<-ggplot(data,aes(x=number,y=V2,group=method))+
    theme_prism()+geom_point(size=2,aes(color=method))+geom_smooth(se=F,aes(color=method),size=1.5)+ylab(i)+xlab("Model evaluation index")+
    scale_color_manual(values = color[1:3])
  ggsave(p,filename = paste0("F:/result/result_11_15/plot/fig3/",str_split(i,":")[[1]][1],"_11_23.pdf"))
}

a<-unique(d[,1])
data<-d[d[,1]==a[3],]
a<-unique(data[,4])
n<-c()
for(i in a){
  data1<-data[data[,4]==i,]
  for(j in 1:9){
    m<-data1[j+1,2]-data1[j,2]
    n<-c(n,m)
  }
}
n<-data.frame(n,c(rep(a[1],9),rep(a[2],9),rep(a[3],9)),rep(c("10-5","15-10","20-15","25-20","30-25","35-30","40-35","45-40","50-45"),3))
colnames(n)<-c("F1","method","number")
ggplot(n,aes(x=number,y=F1,group=method))+
  theme_prism()+geom_point(size=2,aes(color=method))+geom_smooth(se=F,aes(color=method),size=1.5)+ylab("precision")+xlab(" Methylation site interval ")+
  scale_color_manual(values = color[1:3])+geom_hline(aes(yintercept=0.005),color="red")+ theme(axis.text.x = element_text(angle=30, hjust=1, vjust=.9))

a<-list.files("f:/result/result_11_15/process_data/DNN_different_number_11_23/")
d<-data.frame()
for(i in c("DNN_Extra","DNN_Xgboost","DNN_RF")){
  for(j in c(5,10,15,20,25,30,35,40,45,50)){
    b<-read.table(paste0("f:/result/result_11_15/process_data/DNN_different_number_11_23/",i,j,"_11_23.txt"))
    b["number"]<-j
    b["method"]<-i
    #b<-b[5:8,]
    d<-rbind(d,b)
  }
}

a<-unique(d[,1])
for(i in a){
  data<-d[d[,1]==i,]
  p<-ggplot(data,aes(x=number,y=V2,group=method))+
    theme_prism()+geom_point(size=2,aes(color=method))+geom_smooth(se=F,aes(color=method),size=1.5)+ylab(i)+xlab("Model evaluation index")+
    scale_color_manual(values = color[1:3])
  ggsave(p,filename = paste0("F:/result/result_11_15/plot/fig4/",str_split(i,":")[[1]][1],"_11_23.pdf"))
}

a<-unique(d[,1])
data<-d[d[,1]==a[3],]
a<-unique(data[,4])
n<-c()
for(i in a){
  data1<-data[data[,4]==i,]
  for(j in 1:9){
    m<-data1[j+1,2]-data1[j,2]
    n<-c(n,m)
  }
}
n<-data.frame(n,c(rep(a[1],9),rep(a[2],9),rep(a[3],9)),rep(c("10-5","15-10","20-15","25-20","30-25","35-30","40-35","45-40","50-45"),3))
colnames(n)<-c("F1","method","number")
ggplot(n,aes(x=number,y=F1,group=method))+
  theme_prism()+geom_point(size=2,aes(color=method))+geom_smooth(se=F,aes(color=method),size=1.5)+ylab("precision")+xlab(" Methylation site interval ")+
  scale_color_manual(values = color[1:3])+geom_hline(aes(yintercept=0.0),color="red")+ theme(axis.text.x = element_text(angle=30, hjust=1, vjust=.9))

##############################   雷达图
library(pheatmap)
a<-list.files("F:/result/result_11_15/process_data/TCGA_validata_result_11_24//")
b<-grep("_11_24.txt$",a)
a<-a[b]

b<-str_split_fixed(a,"_11_24",2)[,1]

c<-read.table(paste0("F:/result/result_11_15/process_data/TCGA_validata_result_11_24/",a[1]),fill = T)
c<-c[1:4,]
c<-as.data.frame(c)
d<-c

for(i in 2:6){
  c<-read.table(paste0("F:/result/result_11_15/process_data/TCGA_validata_result_11_24/",a[i]),fill = T)
  c<-c[1:4,]
  d<-cbind(d,c)
  }

#d<-d[-c(3,5,7,9,11),]

d<-d[,-c(3,5,7,9,11)]

d<-t(d)
d<-as.data.frame(d)
d["max"]<-1
d["min"]<-0
d<-d[,c(5,6,1:4)]
d<-t(d)
d<-as.data.frame(d)

rownames(d)[3:6]<-d[3:6,1]
d<-d[,-1]
colnames(d)<-b

c<-apply(d, 2, as.numeric)
rownames(c)<-rownames(d)

library(fmsb)
student1_data <- c[c("max", "min", "validata_precision:"), ]
student1_data<-as.data.frame(student1_data)
radarchart(student1_data)

radarchart(
  student1_data, axistype = 1,
  # Customize the polygon
  pcol = "#E7B800", pfcol = scales::alpha("#E7B800", 0.5), plwd = 2, plty = 1,
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 2,
  # Customize the axis
  axislabcol = "grey", 
  # Variable labels
  vlcex = 0.7, vlabels = colnames(student1_data),
  caxislabels = c( 0.5, 0.65,0.75,0.85, 1.0))


