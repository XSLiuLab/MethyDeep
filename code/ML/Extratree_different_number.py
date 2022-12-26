import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score
import sys
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score, matthews_corrcoef,cohen_kappa_score



data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
train_y = data['project_id'].apply(lambda x: Class_dict[x])
train_x = data.iloc[:,1:data.shape[1]]

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_6.csv")
test_y = data['project_id'].apply(lambda x: Class_dict[x])
test_x = data.iloc[:,1:data.shape[1]]


m = sys.argv[1]
print(m)
m=int(m)

Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/Extratree_feature"+str(m)+".csv")

i = sys.argv[2]
print(i)
i=int(i)

from sklearn.ensemble import ExtraTreesClassifier

new_train_x = train_x[Extra_feature.iloc[0:i,0]]
new_test_x = test_x[Extra_feature.iloc[0:i,0]]
clf_Extra = ExtraTreesClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
# scores3 = cross_val_score(clf3, X, y)
# print(scores3.mean())
clf_Extra.fit(new_train_x,train_y)
e_predict = clf_Extra.predict(new_test_x)
accuracy = accuracy_score(test_y,e_predict)
F1_score = f1_score(test_y, e_predict,average='macro')
recall = recall_score(test_y, e_predict, average='weighted')
precision = precision_score(test_y, e_predict,average='weighted')
matthews = matthews_corrcoef(test_y, e_predict,  sample_weight=None)
kappa = cohen_kappa_score(test_y, e_predict)
f = open("/public/slst/home/ningwei/methylation/process_data/ML_different_number_12_15/"+"ExtraTree_"+str(m)+"_"+str(i)+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) 
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) 
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) 
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) 
f.write("\n")   #换行
f.write("test_matthews: "+str(matthews)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_kappa: "+str(kappa)) #将字符串写入文件中
f.close()