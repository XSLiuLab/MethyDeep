from tensorflow.keras.models import load_model
from tensorflow.keras import regularizers
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
import matplotlib.pyplot as plt
from tensorflow.keras.backend import clear_session
from sklearn.model_selection import PredefinedSplit
import math
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import sys
from tensorflow import keras as K
import tensorflow as tf
from tensorflow.keras import regularizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.backend import clear_session
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score, matthews_corrcoef,cohen_kappa_score
import re
from sklearn.ensemble import RandomForestClassifier
###################################
def cal_metrics(confusion_matrix):
  n_classes = confusion_matrix.shape[0]
ppv1 = []
npv1 = []
sen1 = []
spe1 = []
for i in range(n_classes):
  # 逐步获取 真阳，假阳，真阴，假阴四个指标，并计算三个参数
  ALL = np.sum(confusion_matrix)
# 对角线上是正确预测的
TP = confusion_matrix[i, i]
# 列加和减去正确预测是该类的假阳
FP = np.sum(confusion_matrix[:, i]) - TP
# 行加和减去正确预测是该类的假阴
FN = np.sum(confusion_matrix[i, :]) - TP
# 全部减去前面三个就是真阴
TN = ALL - TP - FP - FN
ppv1.append(TP/(TP+FP))
npv1.append(TN/(TN+FN))
sen1.append(TP/(TP+FN))
spe1.append(TN/(TN+FP))
return ppv1, npv1, sen1, spe1

#############################  validata

primary_data = pd.read_csv("/public/slst/home/ningwei/methylation/data/other_model/primary.csv")
primary_type = pd.read_csv("/public/slst/home/ningwei/methylation/data/other_model/primary_cancertype.csv")

metastatic_data = pd.read_csv("/public/slst/home/ningwei/methylation/data/other_model/metastatic.csv")
metastatic_type = pd.read_csv("/public/slst/home/ningwei/methylation/data/other_model/metastatic_cancertype.csv")

#######################
# validata_data = pd.read_csv("/public/slst/home/ningwei/methylation/data/cup_data/test_data.csv")
# validata_really = ["TCGA-HNSC"]*49

data = pd.read_csv("/public/slst/home/ningwei/methylation/data/other_model/train_data.csv")

data['project_id'] = data['project_id'].replace('TCGA-LGG', 'TCGA-GBM') ### 合并LGG_GBM
data['project_id'] = data['project_id'].replace('TCGA-READ', 'TCGA-COAD') ### 合并 READ COAD

data = data.loc[data['project_id']!="TCGA-STAD",]  #### 删除 STAD
data = data.loc[data['project_id']!="TCGA-ESCA",]  #### 删除 ESCA

Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
x = data.iloc[:,1:353]
y = data.iloc[:,0].apply(lambda x: Class_dict[x])

train_x, test_x, train_y, test_y = train_test_split(x, y, stratify=y,\
                                                    train_size=0.8, test_size=0.2, random_state=1)
target_names = list(Class_dict.keys())
####################################### other method
a = os.listdir("/public/slst/home/ningwei/methylation/data/cup_data/")
a = a[2:7]
b = []
for i in range(0,5):
  b.append(a[i].split(".")[0])

for i in range(0,5):
  feature = pd.read_csv("/public/slst/home/ningwei/methylation/data/cup_data/"+str(a[i]))
  new_train_x = train_x[feature.iloc[:,0]]
  clf2 = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
  clf2.fit(new_train_x,train_y)
  #new_validata = metastatic_data[feature.iloc[:,0]]
  new_validata = primary_data[feature.iloc[:,0]]
  R_predict = clf2.predict_proba(new_validata)
  R_predict = pd.DataFrame(R_predict)
  R_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/other_model_result_12_16/primary_ROC/"+str(a[i]),index=0)
  R_predict = clf2.predict(new_validata)
  predict2 = []
  for j in R_predict:
      predict2.append(target_names[j])
  accuracy = accuracy_score(primary_type,predict2)
  F1_score = f1_score(primary_type,predict2,average='weighted')
  recall = recall_score(primary_type, predict2, average='weighted')
  precision = precision_score(primary_type, predict2,average='weighted')
  matthews = matthews_corrcoef(primary_type, predict2,  sample_weight=None)
  kappa = cohen_kappa_score(primary_type, predict2)
  f = open("/public/slst/home/ningwei/methylation/process_data/other_model_result_12_16/primary_"+str(b[i])+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
  f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_recall: "+str(recall)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_precision: "+str(precision)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
  f.write("\n")   #换行
  predict2 = pd.DataFrame(predict2)
  confusion_mat = confusion_matrix(primary_type,predict2,labels=target_names)
  ppv, npv, sen, spe = cal_metrics(confusion_mat)
  f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_npv: "+str(npv)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_sen: "+str(sen)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_spe: "+str(spe)) #将字符串写入文件中
  f.write("\n")   #换行
  f.close()
####################################### DNN
# data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_11_23.csv")
# Class = data['project_id'].unique()
# Class_dict = dict(zip(Class, range(len(Class))))
# train_x = data.iloc[:,1:5078]
# data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(train_y)
train_y = transformed_labels

#xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")
RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/RF_feature45.csv")
target_names = list(Class_dict.keys())

new_train_x = train_x[RF_feature.iloc[0:30,0]]

clear_session()
model = K.models.Sequential()
model.add(K.layers.Dense(units=230, input_dim=30,  activation='relu'))
model.add(K.layers.Dropout(0.1))
model.add(K.layers.Dense(units=192, activation='relu'))
model.add(K.layers.Dropout(0.1))
model.add(K.layers.Dense(units=154, activation='relu'))
model.add(K.layers.Dense(units=116, activation='relu'))
model.add(K.layers.Dense(units=78, activation='relu'))
model.add(K.layers.Dense(units=40, activation='relu'))
model.add(K.layers.Dense(units=26,  activation='softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
b_size = 100
max_epochs = 100

h = model.fit(new_train_x, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
################### mestasis
#new_validata = metastatic_data[RF_feature.iloc[0:30,0]]
new_validata = primary_data[RF_feature.iloc[0:30,0]]
X_predict = model.predict(new_validata)
X_predict = pd.DataFrame(X_predict)
X_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/other_model_result_12_16/primary_ROC/DNN_ROC.csv",index=0)

X_predict = np.argmax(X_predict, axis=-1)
predict2 = []
for j in X_predict: 
  predict2.append(target_names[j])
  accuracy = accuracy_score(primary_type,predict2)
  F1_score = f1_score(primary_type,predict2,average='weighted')
  recall = recall_score(primary_type, predict2, average='weighted')
  precision = precision_score(primary_type, predict2,average='weighted')
  matthews = matthews_corrcoef(primary_type, predict2,  sample_weight=None)
  kappa = cohen_kappa_score(primary_type, predict2)
  f = open("/public/slst/home/ningwei/methylation/process_data/other_model_result_12_16/primary_"+"DNN_RF"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
  f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
  f.write("\n")   #换行 
  f.write("validata_recall: "+str(recall)) #将字符串写入文件中
  f.write("\n")   #换行 
  f.write("validata_precision: "+str(precision)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
  f.write("\n")   #换行
  predict2 = pd.DataFrame(predict2)
  confusion_mat = confusion_matrix(primary_type,predict2,labels=target_names)
  ppv, npv, sen, spe = cal_metrics(confusion_mat)
  f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_npv: "+str(npv)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_sen: "+str(sen)) #将字符串写入文件中
  f.write("\n")   #换行
  f.write("validata_spe: "+str(spe)) #将字符串写入文件中
  f.write("\n")   #换行
  f.close()
