library(readxl)
##############################  PMID: 32415265
a<-read_xlsx("f:/result/result_11_15/other_model/NIHMS1765238-supplement-Supple_table.xlsx",sheet = 10)
a<-as.data.frame(a)
colnames(a)<-1:4
a<-a[-1,]
cg_28<-a[,1]
write.csv(cg_28,file = "F:/result/result_11_15/other_model/cg_28.csv",row.names = F)

a<-read_xlsx("f:/result/result_11_15/other_model/NIHMS1765238-supplement-Supple_table.xlsx",sheet = 11)
a<-as.data.frame(a)
colnames(a)<-1:4
a<-a[-1,]
cg_283<-a[,1]
write.csv(cg_283,file = "F:/result/result_11_15/other_model/cg_283.csv",row.names = F)


a<-read_xlsx("f:/result/result_11_15/other_model/NIHMS1765238-supplement-Supple_table.xlsx",sheet = 12)
a<-as.data.frame(a)
colnames(a)<-1:4
a<-a[-1,]
cg_53<-a[,1]
write.csv(cg_53,file = "F:/result/result_11_15/other_model/cg_53.csv",row.names = F)

######################   PMCID: PMC8770539 

cg_6<-c("cg03292149","cg02573468","cg22427313","cg11056055","cg25978327","cg27598340")
write.csv(cg_6,file = "F:/result/result_11_15/other_model/cg_6.csv",row.names = F)

#################### PMID: 31590287

cg_12<-c("cg01397449","cg04374393","cg06575035","cg07333191","cg16389386","cg16508627",
         "cg16926102","cg17804348","cg19710323","cg22620090","cg26642667","cg26733975")
write.csv(cg_12,file = "F:/result/result_11_15/other_model/cg_12.csv",row.names = F)

################## methdeep

cg_30<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Xgboost_feature_11_23.csv")
cg_30<-cg_30[1:30,1]

all_cg<-c(cg_12,cg_28,cg_283,cg_30,cg_53,cg_6)
all_cg<-unique(all_cg)
write.csv(all_cg,file = "F:/result/result_11_15/other_model/all_cg.csv",row.names = F)

data<-read.csv("F:/result/result_11_15/process_data/test_data.csv")
data<-data[,-1]
a<-str_split_fixed(colnames(data),"\\.",3)[,2]
colnames(data)<-a

b<-setdiff(all_cg,colnames(data))
for(i in b){
  data[i]<-0
}
write.csv(data,file = "F:/result/result_11_15/process_data/test_data.csv",row.names = F)
################### 韦恩图
library(ggvenn)
library(ggVennDiagram)

data<-list(`cg_12` = cg_12,`cg_28`= cg_28,`cg_283` = cg_283,`cg_30` = cg_30,
           `cg_53` = cg_53,`cg_6` = cg_6)

upset(fromList(data),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2, # 文字标签的大小
      main.bar.color=color[1:9],
      sets.bar.color=color[10:15]
      
)

################# 雷达图

a<-list.files("F:/result/result_11_15/process_data/CUP_result/")
a<-a[-7]
b<-str_split_fixed(a,"\\.",2)[,1]

c<-read.table(paste0("F:/result/result_11_15/process_data/CUP_result/",a[1]),fill=T)
c<-c[1:6,]
d<-c

for(i in 2:6){
  c<-read.table(paste0("F:/result/result_11_15/process_data/CUP_result/",a[i]),fill=T)
  c<-c[1:6,]
  d<-cbind(d,c)
  }
d<-d[,-c(3,5,7,9,11)]

d<-t(d)
d<-as.data.frame(d)
d["max"]<-1
d["min"]<-0
d<-d[,c(7,8,1:6)]
d<-t(d)
rownames(d)[3:8]<-d[3:8,1]
d<-d[,-1]
colnames(d)<-b

c<-apply(d, 2, as.numeric)
rownames(c)<-rownames(d)

student1_data <- c[c("max", "min", "validata_acc:"), ]
student1_data<-as.data.frame(student1_data)

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

###########################  生物意义

XGBOOST<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Xgboost_feature_11_23.csv")
Xgboost<-XGBOOST[1:30,1]

data<-fread("F:/result/result_11_15/process_data/illuminaMethyl450_hg38_GDC.csv",data.table = F)
data<-data[data[,1]%in%Xgboost,]

cg_gene<-data[,2]
write.table(cg_gene,file = "F:/result/result_11_15/process_data/cg_gene.txt",quote = F,row.names = F)

df<-read.csv("F:/result/result_11_15/process_data/go_富集/gProfiler_hsapiens_2022-11-25 09-47-04__intersections.csv")
df<-df[1:9,]
df["ratio"]<-df$intersection_size/df$term_size
ggplot(df,aes(x = ratio, 
              y = reorder(term_name,negative_log10_of_adjusted_p_value,sum), # 按照富集度大小排序
              size = intersection_size,
              colour=negative_log10_of_adjusted_p_value)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "Ratio", y = "BP")+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="Enrichment",                        # 图例名称
    low="#3B49927F",                              # 设置颜色范围
    high="#EE00007F")+
  scale_radius(                               # 设置点大小图例
    range=c(4,6),                             # 设置点大小的范围
    name="Size")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()                                  # 设置主题

################### 小提琴图
library(ggplot2)

data<-read.csv("F:/result/result_11_15/process_data/xgboost_data.csv")
a<-colnames(data)[-1]
for(i in 1:30){
  p<-ggplot(data, aes_string(x = "project_id", y = a[i]))+ geom_violin(aes(fill = project_id),trim = FALSE)+theme_classic()+
    scale_fill_manual(values=color)+theme(legend.position = 'none')+ ## 去除legend
    theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))
  ggsave(p,filename = paste0("F:/result/result_11_15/plot/fig6/venn/",a[i],".pdf"))
}

##################### 密度图
XGBOOST<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Xgboost_feature_11_23.csv")
Xgboost<-XGBOOST[1:30,1]
c<-c[c[,1]%in%Xgboost,]

ggplot(c,aes(x=Freq))+geom_density()
ggplot(c,aes(x=Freq))+geom_density(fill="#FC8A61", color=c("#FC8A61"), alpha=0.8)+theme_classic()+
  xlab("number of Recurrent DMS")

################### 附图1
data<-data.frame(c(14174,411,310),c("primary","metastasis","Recurren"))
colnames(data)<-c("number","cancertypes")
data$cancertypes<-factor(data$cancertypes,levels = c("primary","metastasis","Recurren"))
ggplot(data ,aes(x=cancertypes,y=number,fill=cancertypes))+geom_bar(stat = 'identity')+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10)[4:6])+theme_prism()+
  xlab("")+ylab("number of DMS")+theme(legend.position = 'none')+
  geom_text(aes(label=number, y=number+0.05), position=position_dodge(0.9), vjust=0)
#################### 混淆矩阵
##################################   验证集
library(pheatmap)


a<-read.csv("F:/result/result_11_15/process_data/TCGA_validata_result_11_24/xgboost_metastasis_confus_11_24.csv")
a<-a[,-1]
cancertype<-c( 'TCGA-GBM', 'TCGA-CESC', 'TCGA-LUSC', 'TCGA-PCPG', 'TCGA-COAD', 'TCGA-KIRP', 'TCGA-BLCA', 'TCGA-THCA', 'TCGA-BRCA', 'TCGA-UCEC', 'TCGA-LUAD', 'TCGA-THYM', 'TCGA-SARC', 'TCGA-PRAD', 'TCGA-HNSC', 'TCGA-UVM', 'TCGA-PAAD', 'TCGA-MESO', 'TCGA-LAML', 'TCGA-ACC', 'TCGA-LIHC', 'TCGA-KIRC', 'TCGA-SKCM', 'TCGA-KICH', 'TCGA-TGCT', 'TCGA-UCS')
colnames(a)<-rownames(a)<-cancertype

m<-as.data.frame(apply(a, 1, sum))
n<-paste0(rownames(m),"(",m[,1],")")

rownames(a)<-n
cancertype_col<-cancertype
colnames(a)<-cancertype_col
cancertype<-n

############################################################ 计算recall
b<-matrix(ncol = ncol(a), nrow = nrow(a))
for(i in 1:nrow(a)){
  b[i, ]<-as.numeric(a[i, ]/sum(a[i, ]))
}

b1<-b*100
b1<-round(b1)
############################################################other
colnames(b)<-colnames(a)
rownames(b)<-rownames(a)

b2<-b[order(diag(b),decreasing = T),]
b2<-b2[,order(c(diag(b)),decreasing = T)] ######## other data


b3<-b1[order(diag(b),decreasing = T),]
b3<-b3[,order(c(diag(b)),decreasing = T)] ################# 根据recall排序
rownames(b3)<-rownames(b2)
colnames(b3)<-colnames(b2)

# cancertype1<-cancertype_col[order(c(diag(b)),decreasing = T)] ########### 行名和列名
# cancertype2<-cancertype[order(diag(b),decreasing = T)] 
# 
# rownames(b2)<-cancertype2
# colnames(b2)<-cancertype1

pheatmap(b2,cluster_row = FALSE,cluster_col = FALSE,display_numbers = b3,number_format = "%.2f",fontsize=15)

############################################# 计算 precision
c<-matrix(ncol = ncol(a), nrow = nrow(a))
for(i in 1:ncol(a)){
  c[ ,i]<-a[ ,i]/sum(a[ ,i])
}

precis<-diag(c)
precis<-round(100*precis)

annotation_col = data.frame( precision = c(precis))##testdata

rownames(annotation_col) = colnames(a)

annotation_row = data.frame( recall = diag(b))
rownames(annotation_row) = rownames(a)

pheatmap(b2 , annotation_col = annotation_col,  annotation_row = annotation_row, 
         cluster_row = FALSE,  cluster_col = FALSE, display_numbers = b3, number_format = "%.2f", 
         annotation_names_row = TRUE,  annotation_names_col = TRUE,fontsize=15)
#####################
################################################################## 去掉空白癌症

a<-read.csv("F:/result/result_11_15/process_data/TCGA_validata_result_11_24/xgboost_metastasis_confus_11_24.csv")
a<-a[,-1]
#cancertype<-c('TCGA-GBM', 'TCGA-CESC', 'TCGA-LUSC', 'TCGA-PCPG', 'TCGA-COAD', 'TCGA-KIRP', 'TCGA-BLCA', 'TCGA-THCA', 'TCGA-BRCA', 'TCGA-UCEC', 'TCGA-LUAD', 'TCGA-THYM', 'TCGA-SARC', 'TCGA-PRAD', 'TCGA-HNSC', 'TCGA-UVM', 'TCGA-PAAD', 'TCGA-MESO', 'TCGA-LAML', 'TCGA-ACC', 'TCGA-LIHC', 'TCGA-KIRC', 'TCGA-SKCM', 'TCGA-KICH', 'TCGA-TGCT', 'TCGA-UCS')
rownames(a)<-cancertype
colnames(a)<-cancertype_col

m<-which(apply(a, 1,function(x){sum(as.numeric(x))})==0)
m_cancer<-rownames(a)[m]
n<-which(apply(a,2,function(x){sum(as.numeric(x))})==0)
n_cancer<-colnames(a)[n]
#colnames(b3)<-rownames(b3)<-colnames(b2)


b2<-b2[!rownames(b2)%in%m_cancer,!colnames(b2)%in%n_cancer]

# cancer_number<-read.csv("F:/result/validata_plot/cancer_number.csv")
# for(i in 1:nrow(b2)){
#   m<-which(cancer_number[,1]==rownames(b2)[i])
#   if(length(m)==1){
#     rownames(b2)[i]<-paste0(rownames(b2)[i],"_",cancer_number[m,2],")")
#   }
# }
# rownames(b2)<-str_split_fixed(rownames(b2),"\\(",2)[,1]
pheatmap(b2 , annotation_col = annotation_col[!annotation_col%in%n_cancer],  annotation_row = annotation_row[!annotation_row%in%m_cancer], 
         cluster_row = FALSE,  cluster_col = FALSE, display_numbers = b3[!rownames(b3)%in%m_cancer,!colnames(b3)%in%n_cancer], number_format = "%.2f", 
         annotation_names_row = TRUE,  annotation_names_col = TRUE,fontsize=15)
####保存PDF为 high 9.47 width 14.54 inches

######################## 条形图
m<-list.files("F:/result/result_11_15/process_data/CUP_result/")
m<-m[-7]
n<-str_split_fixed(m,"\\.",2)[,1]  
d<-data.frame()
for(j in 1:length(m)){
  data<-read.table(paste0("F:/result/result_11_15/process_data/CUP_result/",m[j]),fill = T,sep = "\t")
  data<-data[7:10,]
  index<-c("PPV","NPV","SEN","SPE")
  c<-data.frame()
  for(i in 1:4){
    a<-str_split_fixed(data[i],"\\,",26)
    a[1]<-str_split_fixed(a[1],"\\[",2)[,2]
    a[26]<-str_split_fixed(a[26],"\\]",2)[,1]
    b<-data.frame(t(a),rep(index[i],26))
    c<-rbind(c,b)
  }
  c[,3]<-n[j]
  d<-rbind(c,d)
}
d["cancer"]<-rep(cancertype_col,24)
d[,1]<-as.numeric(d[,1])
colnames(d)[1:3]<-c("number","index","method")

d1<-d[d[,4]=="TCGA-COAD",]

d1[,3]<-factor(d1[,3],levels = c("DNN_xgboost","cg_283","cg_53","cg_28","cg_12","cg_6"))
ggplot(d1,aes(x=method,y=number,fill=index))+geom_bar(stat = 'identity',position = position_dodge(0.9))+
  scale_fill_manual( values = pal_npg("nrc", alpha =0.7)(10))+theme_prism()+theme(axis.text.x = element_text(angle=30, hjust=0.9, vjust=.9))+
  xlab("")+ylab("number of DMS")

+theme(legend.position = 'none') ## 去除legend
