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
########################################################################### 构建划分函数

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

#train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote_11_22.csv")

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_11_23.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
# train_y = []
# for i in range(len(data.iloc[:,0])):
#     if data.iloc[i,0] in Class_dict.keys():
#         train_y.append(Class_dict[data.iloc[i,0]])
train_x = data.iloc[:,1:5078]
data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(data['target'])
train_y = transformed_labels

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_11_23.csv")
# test_y = []
# for i in range(len(data.iloc[:,0])):
#     if data.iloc[i,0] in Class_dict.keys():
#         test_y.append(Class_dict[data.iloc[i,0]])
test_x = data.iloc[:,1:5078]
data["target"] = data.iloc[:,0].apply(lambda x: Class_dict[x])
# 对目标变量进行0-1编码(One-hot Encoding)
lb = LabelBinarizer()
lb.fit(list(Class_dict.values()))
transformed_labels = lb.transform(data['target'])
test_y = transformed_labels

RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/RF_feature_11_23.csv")
Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Extratree_feature_11_23.csv")
xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")

m = sys.argv[1]
m =int(m)

########################## Extra
new_train_x = train_x[Extra_feature.iloc[0:m,0]]
new_test_x = test_x[Extra_feature.iloc[0:m,0]]
if m ==5 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=50, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=48, activation='relu'))
    model.add(K.layers.Dense(units=46, activation='relu'))
    model.add(K.layers.Dense(units=44, activation='relu'))
    model.add(K.layers.Dense(units=42, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100
if m ==10 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=100, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=88, activation='relu'))
    model.add(K.layers.Dense(units=75, activation='relu'))
    model.add(K.layers.Dense(units=63, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==15 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=130, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=114, activation='relu'))
    model.add(K.layers.Dense(units=98, activation='relu'))
    model.add(K.layers.Dense(units=82, activation='relu'))
    model.add(K.layers.Dense(units=66, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==20 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=200, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=168, activation='relu'))
    model.add(K.layers.Dense(units=136, activation='relu'))
    model.add(K.layers.Dense(units=104, activation='relu'))
    model.add(K.layers.Dense(units=72, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==25 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=250, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=198, activation='relu'))
    model.add(K.layers.Dense(units=145, activation='relu'))
    model.add(K.layers.Dense(units=93, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100
if m ==30 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=250, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=200, activation='relu'))
    model.add(K.layers.Dense(units=150, activation='relu'))
    model.add(K.layers.Dense(units=100, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==35 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=310, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=243, activation='relu'))
    model.add(K.layers.Dense(units=175, activation='relu'))
    model.add(K.layers.Dense(units=108, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m ==40 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=340, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=268, activation='relu'))
    model.add(K.layers.Dense(units=195, activation='relu'))
    model.add(K.layers.Dense(units=123, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 100
    max_epochs = 100
if m ==45 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=400, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=313, activation='relu'))
    model.add(K.layers.Dense(units=225, activation='relu'))
    model.add(K.layers.Dense(units=138, activation='relu'))
    model.add(K.layers.Dense(units=50, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100
if m == 50 :
    clear_session()
    model = K.models.Sequential()
    model.add(K.layers.Dense(units=430, input_dim=m,  activation='relu'))
    model.add(K.layers.Dense(units=300, activation='relu'))
    model.add(K.layers.Dense(units=170, activation='relu'))
    model.add(K.layers.Dense(units=40, activation='relu'))
    model.add(K.layers.Dense(units=26,  activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    b_size = 80
    max_epochs = 100

h = model.fit(new_train_x, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数

R_predict = model.predict(new_test_x)
R_predict = np.argmax(R_predict, axis=-1)

y_really = np.argmax(test_y, axis=-1)

accuracy = accuracy_score(y_really,R_predict)
F1_score = f1_score(y_really,R_predict,average='weighted')
recall = recall_score(y_really, R_predict, average='weighted')
precision = precision_score(y_really, R_predict,average='weighted')
f = open("/public/slst/home/ningwei/methylation/process_data/DNN_different_number/"+"DNN_Extra"+str(m)+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
f.write("test_acc: "+str(accuracy)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_recall: "+str(recall)) #将字符串写入文件中
f.write("\n")   #换行 
f.write("test_precision: "+str(precision)) #将字符串写入文件中
f.write("\n")   #换行
f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
f.write("\n")   #换行
f.close()