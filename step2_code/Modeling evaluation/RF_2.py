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
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score

gred_list = [d for d in os.listdir('/public/slst/home/ningwei/TCGA/Methylation450K/data_6_15')]
data = pd.read_csv("~/TCGA/Methylation450K/data_6_15/"+gred_list[2])
target_names = data.iloc[:,0].unique()

for m in [5,10,15,20,25,30,35,40,45,50]:
    clear_session()
    train_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_x_2"+"_"+str(m)+".csv")
    train_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_y_2"+"_"+str(m)+".csv")
    test_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_test_x_2"+"_"+str(m)+".csv")
    test_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_test_y_2"+"_"+str(m)+".csv")
    ###############################################################  对数据进行标准化
    train_x_zcore=(train_x-train_x.mean(axis=0))/train_x.std(axis=0)
    test_x_zcore=(test_x-train_x.mean(axis=0))/train_x.std(axis=0)
    if m==5:
        b_size = 80
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=30, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dense(units=30, activation='relu'))
        model.add(K.layers.Dense(units=30, activation='relu'))
        model.add(K.layers.Dense(units=30, activation='relu'))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==10:
        b_size = 80
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=100, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.3))
        model.add(K.layers.Dense(units=63, activation='relu'))
        model.add(K.layers.Dropout(0.15))
        model.add(K.layers.Dense(units=25, activation='relu'))
        model.add(K.layers.Dropout(0))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==15:
        b_size = 80
        max_epochs = 50
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=150, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dense(units=109, activation='relu'))
        model.add(K.layers.Dense(units=67, activation='relu'))
        model.add(K.layers.Dense(units=26, activation='relu'))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==20:
        b_size = 50
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=150, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.3))
        model.add(K.layers.Dense(units=109, activation='relu'))
        model.add(K.layers.Dropout(0.2))
        model.add(K.layers.Dense(units=67, activation='relu'))
        model.add(K.layers.Dropout(0.1))
        model.add(K.layers.Dense(units=26, activation='relu'))
        model.add(K.layers.Dropout(0.0))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==25:
        b_size = 100
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=230, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.5))
        model.add(K.layers.Dense(units=180, activation='relu'))
        model.add(K.layers.Dropout(0.38))
        model.add(K.layers.Dense(units=130, activation='relu'))
        model.add(K.layers.Dropout(0.26))
        model.add(K.layers.Dense(units=80, activation='relu'))
        model.add(K.layers.Dropout(0.14))
        model.add(K.layers.Dense(units=30, activation='relu'))
        model.add(K.layers.Dropout(0.02))  
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==30:
        b_size = 50
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=300, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.3))
        model.add(K.layers.Dense(units=165, activation='relu'))
        model.add(K.layers.Dropout(0.15))
        model.add(K.layers.Dense(units=30, activation='relu'))
        model.add(K.layers.Dropout(0.0))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==35:
        b_size = 50
        max_epochs = 50
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=250, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.3))
        model.add(K.layers.Dense(units=194, activation='relu'))
        model.add(K.layers.Dropout(0.23))
        model.add(K.layers.Dense(units=138, activation='relu'))
        model.add(K.layers.Dropout(0.16))
        model.add(K.layers.Dense(units=82, activation='relu'))
        model.add(K.layers.Dropout(0.09))
        model.add(K.layers.Dense(units=25, activation='relu'))
        model.add(K.layers.Dropout(0.02))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==40:
        b_size = 80
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=280, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.5))
        model.add(K.layers.Dense(units=153, activation='relu'))
        model.add(K.layers.Dropout(0.25))
        model.add(K.layers.Dense(units=25, activation='relu'))
        model.add(K.layers.Dropout(0.0))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==45:
        b_size = 100
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=410, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.5))
        model.add(K.layers.Dense(units=314, activation='relu'))
        model.add(K.layers.Dropout(0.23))
        model.add(K.layers.Dense(units=218, activation='relu'))
        model.add(K.layers.Dropout(0.16))
        model.add(K.layers.Dense(units=122, activation='relu'))
        model.add(K.layers.Dropout(0.09))
        model.add(K.layers.Dense(units=25, activation='relu'))
        model.add(K.layers.Dropout(0.02))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    if m==50:
        b_size = 50
        max_epochs = 100
        clear_session()
        model = K.models.Sequential()
        model.add(K.layers.Dense(units=280, input_dim=train_x_zcore.shape[1],  activation='relu'))
        model.add(K.layers.Dropout(0.5))
        model.add(K.layers.Dense(units=153, activation='relu'))
        model.add(K.layers.Dropout(0.25))
        model.add(K.layers.Dense(units=25, activation='relu'))
        model.add(K.layers.Dropout(0.00))
        model.add(K.layers.Dense(units=21,  activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    h = model.fit(train_x_zcore, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
    actual = np.argmax(test_y.values, axis=-1)
    predicted = np.argmax(model.predict(test_x_zcore), axis=-1)
    acc = accuracy_score(actual, predicted)
    recall = recall_score(actual, predicted, average='weighted')
    precision = precision_score(actual,predicted,average='weighted')
    f = open("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/finall_result/"+"new_2"+"_"+str(m)+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行 
    other_data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/random_2.csv")
    other_data1 = copy.deepcopy(other_data)
    other_data1.columns = other_data1.columns+".1"
    other_data = pd.concat([other_data,other_data1],axis=1)
    other_data = other_data[train_x_zcore.columns]
    other_data_zscore = (other_data-train_x.mean(axis=0))/train_x.std(axis=0)
    predict = np.argmax(model.predict(other_data_zscore), axis=-1)
    predict2 = []
    for h in predict: 
        predict2.append(target_names[h])
    other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/random_2_cancertype.csv")
    acc = accuracy_score(other_really, predict2)
    recall = recall_score(other_really, predict2, average='weighted')
    precision = precision_score(other_really, predict2,average='weighted')
    F1_score=f1_score(other_really,predict2,average="micro")
    f.write("vali_acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("vali_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("vali_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("vali_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()