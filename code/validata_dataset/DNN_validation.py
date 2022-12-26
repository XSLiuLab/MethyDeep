#########################################################  加载模块
import sys
from tensorflow import keras as K
import tensorflow as tf
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
from tensorflow.keras.models import load_model
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
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
########################################################################### 构建划分函数



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

#train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote.csv")

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(data['target'])
train_y = transformed_labels
train_x = data.iloc[:,1:data.shape[1]]

# data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_7.csv")
# data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# # 对目标变量进行0-1编码(One-hot Encoding)
# lb = LabelBinarizer()
# lb.fit(list(Class_dict.values()))
# transformed_labels = lb.transform(data['target'])
# test_y = transformed_labels
# test_x = data.iloc[:,1:data.shape[1]]

RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/RF_feature45.csv")
# Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Extratree_feature_11_23.csv")
# xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")

# validata_metastasis = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata/metastatic.csv")
# n = validata_metastasis.iloc[:,0]!="TCGA-ESCA"
# validata_metastasis = validata_metastasis.loc[n,]
# metastatic_validata_really = validata_metastasis.iloc[:,0]
#validata_Recurren = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata/recurren.csv")
#recurren_validata_really = validata_Recurren.iloc[:,0]

# #############################    RF validata
# geo_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/GEO_metasta_validata_cancertype.csv")
# geo_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/GEO_metasta_validata.csv")
# metastatic_validata_really = pd.DataFrame(metastatic_validata_really)
# geo_metasta_validata_really.columns = metastatic_validata_really.columns
# RF_metastatic_validata_really = pd.concat([metastatic_validata_really,geo_metasta_validata_really])
# #############################    Extratree validata
# EXtra_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/Extratree_GEO_metasta_validata_cancertype.csv")
# Extra_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/Extratree_GEO_metasta_validata.csv")
# EXtra_metasta_validata_really.columns = metastatic_validata_really.columns
# EXtra_metasta_validata_really = pd.concat([metastatic_validata_really,EXtra_metasta_validata_really])
# #############################    xgboost validata
# xgboost_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/xgboost_GEO_metasta_validata_cancertype.csv")
# xgboost_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/xgboost_GEO_metasta_validata.csv")
# xgboost_metasta_validata_really.columns = metastatic_validata_really.columns
# xgboost_metasta_validata_really = pd.concat([metastatic_validata_really,xgboost_metasta_validata_really])

# #data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote.csv")
# #target_names = data.iloc[:,0].unique()
target_names = list(Class_dict.keys())
########################## RF 
new_train_x = train_x[RF_feature.iloc[0:30,0]]
#new_test_x = test_x[RF_feature.iloc[0:25,0]]

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
model.save('/public/slst/home/ningwei/methylation/methydeep_30')
# target_names = pd.DataFrame(target_names)
# target_names.to_csv("/public/slst/home/ningwei/methylation/process_data/targer_names_12_9.csv",index=0)

# #############
# model = load_model('/public/slst/home/ningwei/methylation/methydeep_25')  

# target_names = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/targer_names_12_9.csv")
# target_names = target_names.iloc[:,0]
# validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/validata/primary_oversample.csv")
# new_validata = validata[RF_feature.iloc[0:25,0]]
# validata_type = pd.read_csv("/public/slst/home/ningwei/methylation/data/validata/primary_oversample_cancertype.csv")

# R_predict = model.predict(new_validata)
# R_predict = np.argmax(R_predict, axis=-1)
# predict2 = []
# for j in R_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(validata_type,predict2)
# F1_score = f1_score(validata_type,predict2,average='weighted')
# recall = recall_score(validata_type, predict2, average='weighted')
# precision = precision_score(validata_type, predict2,average='weighted')
# matthews = matthews_corrcoef(validata_type, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(validata_type, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result_12_9/"+"DNN_RF_25"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(validata_type,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result_12_9/DNN_RF_25_confus.csv")

# from sklearn.metrics import classification_report

# print(classification_report(validata_type, predict2))



# ################# RF metastasis
# RF_data = validata_metastasis[RF_feature.iloc[0:30,0]]  ### TCGA-metastasic
# RF_data1 = geo_metasta_validata[RF_feature.iloc[0:30,0]] ### GEO-metastasic
# RF_data = pd.concat([RF_data,RF_data1])
# R_predict = model.predict(RF_data)
# R_predict = pd.DataFrame(R_predict)
# R_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/DNN_RF_metasta.csv",index=0)

# R_predict = np.argmax(R_predict, axis=-1)
# predict2 = []
# for j in R_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(RF_metastatic_validata_really,predict2)
# F1_score = f1_score(RF_metastatic_validata_really,predict2,average='weighted')
# recall = recall_score(RF_metastatic_validata_really, predict2, average='weighted')
# precision = precision_score(RF_metastatic_validata_really, predict2,average='weighted')
# matthews = matthews_corrcoef(RF_metastatic_validata_really, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(RF_metastatic_validata_really, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_RF_metastasis_11_24"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(RF_metastatic_validata_really,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_RF_metastasis_confus_11_24.csv")
# ################### RF recurren
# RF_data = validata_Recurren[RF_feature.iloc[0:30,0]]
# R_predict = model.predict(RF_data)
# R_predict = np.argmax(R_predict, axis=-1)
# predict2 = []
# for j in R_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(recurren_validata_really,predict2)
# F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
# recall = recall_score(recurren_validata_really, predict2, average='weighted')
# precision = precision_score(recurren_validata_really, predict2,average='weighted')
# matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(recurren_validata_really, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_RF_recurren_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_RF_recurren_confus_11_23.csv")
# ########################## Extratree
# new_train_x = train_x[Extra_feature.iloc[0:30,0]]
# new_test_x = test_x[Extra_feature.iloc[0:30,0]]
# clear_session()
# model = K.models.Sequential()
# model.add(K.layers.Dense(units=250, input_dim=30,  activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=200, activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=150, activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=100, activation='relu'))
# model.add(K.layers.Dropout(0.1))
# model.add(K.layers.Dense(units=50, activation='relu'))
# model.add(K.layers.Dense(units=26,  activation='softmax'))
# model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

# b_size = 80
# max_epochs = 100
# h = model.fit(new_train_x, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
# ######################### EXtra metastasis
# Extratree_data = validata_metastasis[Extra_feature.iloc[0:30,0]]
# Extratree_data1 = Extra_metasta_validata[Extra_feature.iloc[0:30,0]] ### GEO-metastasic
# Extratree_data = pd.concat([Extratree_data,Extratree_data1])
# E_predict = model.predict(Extratree_data)
# E_predict = pd.DataFrame(E_predict)
# E_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/DNN_Extra_metasta.csv",index=0)

# E_predict = np.argmax(E_predict, axis=-1)
# predict2 = []
# for j in E_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(EXtra_metasta_validata_really,predict2)
# F1_score = f1_score(EXtra_metasta_validata_really,predict2,average='weighted')
# recall = recall_score(EXtra_metasta_validata_really, predict2, average='weighted')
# precision = precision_score(EXtra_metasta_validata_really, predict2,average='weighted')
# matthews = matthews_corrcoef(EXtra_metasta_validata_really, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(EXtra_metasta_validata_really, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_Extra_metastasis"+"_11_24.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(EXtra_metasta_validata_really,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_Extra_metastasis_confus_11_24.csv")
# # ################### recurren
# Extratree_data = validata_Recurren[Extra_feature.iloc[0:30,0]]
# E_predict = model.predict(Extratree_data)
# E_predict = np.argmax(E_predict, axis=-1)
# predict2 = []
# for j in E_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(recurren_validata_really,predict2)
# F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
# recall = recall_score(recurren_validata_really, predict2, average='weighted')
# precision = precision_score(recurren_validata_really, predict2,average='weighted')
# matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(recurren_validata_really, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_Extra_recurren"+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_Extra_recurren_confus_11_23.csv")
######################### xgboost
# new_train_x = train_x[xgboost_feature.iloc[0:30,0]]
# new_test_x = test_x[xgboost_feature.iloc[0:30,0]]
# clear_session()
# model = K.models.Sequential()
# model.add(K.layers.Dense(units=250, input_dim=30,  activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=210, activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=170, activation='relu'))
# model.add(K.layers.Dropout(0.3))
# model.add(K.layers.Dense(units=130, activation='relu'))
# model.add(K.layers.Dropout(0.1))
# model.add(K.layers.Dense(units=90, activation='relu'))
# model.add(K.layers.Dense(units=50, activation='relu'))
# model.add(K.layers.Dense(units=26,  activation='softmax'))
# model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

# b_size = 100
# max_epochs = 100
# h = model.fit(new_train_x, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数

# model.save('/public/slst/home/ningwei/methylation/methydeep')



# ################### mestasis
# xgboost_data = validata_metastasis[xgboost_feature.iloc[0:30,0]]
# xgboost_data1 = xgboost_metasta_validata[xgboost_feature.iloc[0:30,0]] ### GEO-metastasic
# xgboost_data = pd.concat([xgboost_data,xgboost_data1])
# X_predict = model.predict(xgboost_data)
# X_predict = pd.DataFrame(X_predict)
# X_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/DNN_xgboost_metasta.csv",index=0)


# X_predict = np.argmax(X_predict, axis=-1)
# predict2 = []
# for j in X_predict: 
#     predict2.append(target_names[j])
# accuracy = accuracy_score(xgboost_metasta_validata_really,predict2)
# F1_score = f1_score(xgboost_metasta_validata_really,predict2,average='weighted')
# recall = recall_score(xgboost_metasta_validata_really, predict2, average='weighted')
# precision = precision_score(xgboost_metasta_validata_really, predict2,average='weighted')
# matthews = matthews_corrcoef(xgboost_metasta_validata_really, predict2,  sample_weight=None)
# kappa = cohen_kappa_score(xgboost_metasta_validata_really, predict2)
# f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_xgboost_metastasis"+"_11_24.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# f.write("\n")   #换行 
# f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# f.write("\n")   #换行
# predict2 = pd.DataFrame(predict2)
# confusion_mat = confusion_matrix(xgboost_metasta_validata_really,predict2,labels=target_names)
# ppv, npv, sen, spe = cal_metrics(confusion_mat)
# f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# f.write("\n")   #换行
# f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# f.write("\n")   #换行
# f.close()
# confusion_mat1 = pd.DataFrame(confusion_mat)
# confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_xgboost_metastasis_confus_11_24.csv")
# # ################### recurren
# # xgboost_data = validata_Recurren[xgboost_feature.iloc[0:30,0]]
# # X_predict = model.predict(xgboost_data)
# # X_predict = np.argmax(X_predict, axis=-1)
# # predict2 = []
# # for j in X_predict: 
# #     predict2.append(target_names[j])
# # accuracy = accuracy_score(recurren_validata_really,predict2)
# # F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
# # recall = recall_score(recurren_validata_really, predict2, average='weighted')
# # precision = precision_score(recurren_validata_really, predict2,average='weighted')
# # matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
# # kappa = cohen_kappa_score(recurren_validata_really, predict2)
# # f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"DNN_xgboost_recurren"+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
# # f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
# # f.write("\n")   #换行 
# # f.write("validata_recall: "+str(recall)) #将字符串写入文件中
# # f.write("\n")   #换行 
# # f.write("validata_precision: "+str(precision)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
# # f.write("\n")   #换行
# # predict2 = pd.DataFrame(predict2)
# # confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
# # ppv, npv, sen, spe = cal_metrics(confusion_mat)
# # f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_npv: "+str(npv)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_sen: "+str(sen)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.write("validata_spe: "+str(spe)) #将字符串写入文件中
# # f.write("\n")   #换行
# # f.close()
# # confusion_mat1 = pd.DataFrame(confusion_mat)
# # confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/DNN_xgboost_recurren_confus_11_23.csv")