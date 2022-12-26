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
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score
from sklearn.metrics import matthews_corrcoef,cohen_kappa_score

########################################################################### 构建划分函数

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(data['target'])
train_y = transformed_labels
train_x = data.iloc[:,1:data.shape[1]-1]

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_6.csv")
data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(data['target'])
test_y = transformed_labels
test_x = data.iloc[:,1:data.shape[1]-1]

n = sys.argv[1]
n = int(n)
print(n)
RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/RF_feature"+str(n)+".csv")
# Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Extratree_feature_11_23.csv")
# xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")

m = sys.argv[2]
m =int(m)
print(m)
########################## RF 
new_train_x = train_x[RF_feature.iloc[0:m,0]]
new_test_x = test_x[RF_feature.iloc[0:m,0]]
if m ==5 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=50, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=47, activation='relu'))
    model.add(K.layers.Dropout(0.07))
    model.add(K.layers.Dense(units=44, activation='relu'))
    model.add(K.layers.Dropout(0.04))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dropout(0.01))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==10 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=100, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=90, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=80, activation='relu'))
    model.add(K.layers.Dense(units=70, activation='relu'))
    model.add(K.layers.Dense(units=60, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==15 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=130, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=110, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=90, activation='relu'))
    model.add(K.layers.Dense(units=70, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==20 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=180, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=152, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=124, activation='relu'))
    model.add(K.layers.Dense(units=96, activation='relu'))
    model.add(K.layers.Dense(units=68, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==25 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=250, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=210, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=170, activation='relu'))
    model.add(K.layers.Dense(units=130, activation='relu'))
    model.add(K.layers.Dense(units=90, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==30 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=230, input_dim=m,  activation='relu'))
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
if m ==35 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=330, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=274, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=218, activation='relu'))
    model.add(K.layers.Dense(units=162, activation='relu'))
    model.add(K.layers.Dense(units=106, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100
if m ==40 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=310, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=245, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=180, activation='relu'))
    model.add(K.layers.Dense(units=115, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100
if m ==45 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=400, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=284, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=167, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m == 50 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=310, input_dim=m,  activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=220, activation='relu'))
    model.add(K.layers.Dropout(0.1))
    model.add(K.layers.Dense(units=130, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100

h = model.fit(new_train_x, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数


R_predict = model.predict(new_test_x)
R_predict = np.argmax(R_predict, axis=-1)

y_really = np.argmax(test_y, axis=-1)

accuracy = accuracy_score(y_really,R_predict)
F1_score = f1_score(y_really,R_predict,average='weighted')
recall = recall_score(y_really, R_predict, average='weighted')
precision = precision_score(y_really, R_predict,average='weighted')

matthews = matthews_corrcoef(y_really, R_predict,  sample_weight=None)
kappa = cohen_kappa_score(y_really, R_predict)

f = open("/public/slst/home/ningwei/methylation/process_data/DNN_different_number_12_15/"+"DNN_RF_"+str(n)+"_"+str(m)+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_matthews: "+str(matthews)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_kappa: "+str(kappa)) #将字符串写入文件中
f.close()
