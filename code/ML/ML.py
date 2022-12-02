import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score

#######################################
def load_data(CSV_FILE_PATH):
    IRIS = pd.read_csv(CSV_FILE_PATH)
    target_var = 'project_id'  # 目标变量
    # 数据集的特征
    features = list(IRIS.columns)
    features.remove(target_var)
    # 目标变量的类别
    Class = IRIS[target_var].unique()
    # 目标变量的类别字典
    Class_dict = dict(zip(Class, range(len(Class))))
    # 增加一列target, 将目标变量进行编码
    IRIS['target'] = IRIS[target_var].apply(lambda x: Class_dict[x])
    # 对目标变量进行0-1编码(One-hot Encoding)
    lb = LabelBinarizer()
    lb.fit(list(Class_dict.values()))
    transformed_labels = lb.transform(IRIS['target'])
    y_bin_labels = []  # 对多分类进行0-1编码的变量
    for i in range(transformed_labels.shape[1]):
        y_bin_labels.append('y' + str(i))
        IRIS['y' + str(i)] = transformed_labels[:, i]
    # 将数据集分为训练集和测试集
    train_x, test_x, train_y, test_y = train_test_split(IRIS[features], IRIS[y_bin_labels], stratify=IRIS['project_id'],\
                                                        train_size=0.8, test_size=0.2, random_state=1)
    return train_x, test_x, train_y, test_y

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_11_23.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
train_y = []
for i in range(len(data.iloc[:,0])):
    if data.iloc[i,0] in Class_dict.keys():
        train_y.append(Class_dict[data.iloc[i,0]])
train_x = data.iloc[:,1:5078]
data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_11_23.csv")
test_y = []
for i in range(len(data.iloc[:,0])):
    if data.iloc[i,0] in Class_dict.keys():
        test_y.append(Class_dict[data.iloc[i,0]])
test_x = data.iloc[:,1:5078]

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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
f.write("\n")   #换行
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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
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
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
f.write("\n")   #换行
f.close()
