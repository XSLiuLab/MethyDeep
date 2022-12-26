import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score, matthews_corrcoef,cohen_kappa_score
import sys

#######################################
data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
train_y = data['project_id'].apply(lambda x: Class_dict[x])
train_x = data.iloc[:,1:data.shape[1]]

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_6.csv")
test_y = data['project_id'].apply(lambda x: Class_dict[x])
test_x = data.iloc[:,1:data.shape[1]]

target_names = list(Class_dict.keys())

##########################  K-NN
import pandas as pd
import numpy as np
import operator
import os
from sklearn.neighbors import KNeighborsClassifier

clf_sk = KNeighborsClassifier()
clf_sk.fit(train_x, train_y)



k_predcit = clf_sk.predict(test_x)
accuracy = accuracy_score(test_y,k_predcit)
F1_score = f1_score(test_y, k_predcit,average='macro')
recall = recall_score(test_y, k_predcit, average='weighted')
precision = precision_score(test_y, k_predcit,average='weighted')

f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"K-紧邻_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   #换行
f.close()

################### NB
from sklearn.naive_bayes  import MultinomialNB
nb=MultinomialNB()
nb.fit(train_x, train_y)
print(nb.predict(test_x))

p_predict = nb.predict(test_x)
accuracy = accuracy_score(test_y,p_predict)
F1_score = f1_score(test_y, p_predict,average='macro')
recall = recall_score(test_y, p_predict, average='weighted')
precision = precision_score(test_y, p_predict,average='weighted')

f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"朴素贝叶斯_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   #换行
f.close()

################### 
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier

## DT
clf1 = DecisionTreeClassifier(max_depth=None, min_samples_split=2,random_state=0)
#scores1 = cross_val_score(clf1, train_x, train_y)
#print(scores1.mean())
clf1.fit(train_x,train_y)
print(clf1.predict(test_x))

j_predict = clf1.predict(test_x)
accuracy = accuracy_score(test_y,j_predict)
F1_score = f1_score(test_y, j_predict,average='macro')
recall = recall_score(test_y, j_predict, average='weighted')
precision = precision_score(test_y, j_predict,average='weighted')

f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"决策树_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   #换行
f.close()


## RF
clf2 = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
#scores2 = cross_val_score(clf2, X, y)
#print(scores2.mean())
clf2.fit(train_x,train_y)
print(clf2.predict(test_x))
R_predict = clf2.predict(test_x)
accuracy = accuracy_score(test_y,R_predict)
F1_score = f1_score(test_y,R_predict,average='macro')
recall = recall_score(test_y, R_predict, average='weighted')
precision = precision_score(test_y, R_predict,average='weighted')
f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"随机森林_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   
f.close()

## ExtraTree
clf_Extra = ExtraTreesClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
# scores3 = cross_val_score(clf3, X, y)
# print(scores3.mean())
clf_Extra.fit(train_x,train_y)
print(clf_Extra.predict(test_x))
e_predict = clf_Extra.predict(test_x)
accuracy = accuracy_score(test_y,e_predict)
F1_score = f1_score(test_y, e_predict,average='macro')
recall = recall_score(test_y, e_predict, average='weighted')
precision = precision_score(test_y, e_predict,average='weighted')
f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"ExtraTree_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   #换行
f.close()

###################  XGBOOST
import xgboost as xgb

model = xgb.XGBClassifier(max_depth=5,learning_rate=0.1,n_estimators=160,silent=True,objective='multi:softmax')
model.fit(train_x,train_y)
X_pred = model.predict(test_x)
accuracy = accuracy_score(test_y,X_pred)
F1_score = f1_score(test_y, X_pred,average='macro')
recall = recall_score(test_y, X_pred, average='weighted')
precision = precision_score(test_y, X_pred,average='weighted')
f = open("/public/slst/home/ningwei/methylation/process_data/ML_result/"+"XGBOOST_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) # 
f.write("\n")   # 
f.write("test_recall: "+str(recall)) # 
f.write("\n")   # 
f.write("test_precision: "+str(precision)) # 
f.write("\n")   
f.write("test_F1: "+str(F1_score)) # 
f.write("\n")   
f.close()
