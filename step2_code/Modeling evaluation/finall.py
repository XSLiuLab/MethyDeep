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
import os
import copy

gred_list = [d for d in os.listdir('/public/slst/home/ningwei/TCGA/Methylation450K/data_6_15')]
for j in range(1,100):
    for i in [50]:
        for m in [5,10,15,20,25,30]:
            data = pd.read_csv("~/TCGA/Methylation450K/data_6_15/"+gred_list[i])
            target_names = data.iloc[:,0].unique()
            train_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/new_RF_train_x_"+str(i)+"_"+str(m)+".csv")
            train_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/new_RF_train_y_"+str(i)+"_"+str(m)+".csv")
            test_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/new_RF_test_x_"+str(i)+"_"+str(m)+".csv")
            test_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/new_RF_test_y_"+str(i)+"_"+str(m)+".csv")
            ###############################################################  对数据进行标准化
            train_x_zcore=(train_x-train_x.mean(axis=0))/train_x.std(axis=0)
            test_x_zcore=(test_x-train_x.mean(axis=0))/train_x.std(axis=0)
            
            if m==5:
                b_size = 50
                max_epochs = 100
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=50, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dense(units=45, activation='relu'))      
                model.add(K.layers.Dense(units=40, activation='relu'))
                model.add(K.layers.Dense(units=35, activation='relu'))   
                model.add(K.layers.Dense(units=30, activation='relu'))  
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            if m==10:
                b_size = 80
                max_epochs = 50
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=80, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dense(units=69, activation='relu'))
                model.add(K.layers.Dense(units=58, activation='relu'))
                model.add(K.layers.Dense(units=47, activation='relu'))
                model.add(K.layers.Dense(units=36, activation='relu'))
                model.add(K.layers.Dense(units=25, activation='relu'))
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            if m==15:
                b_size = 80
                max_epochs = 100
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=150, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dense(units=126, activation='relu'))
                model.add(K.layers.Dense(units=102, activation='relu'))
                model.add(K.layers.Dense(units=78, activation='relu'))
                model.add(K.layers.Dense(units=54, activation='relu'))
                model.add(K.layers.Dense(units=30, activation='relu'))
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            if m==20:
                b_size = 100
                max_epochs = 100
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=130, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dense(units=110,  activation='relu'))
                model.add(K.layers.Dense(units=90,  activation='relu'))
                model.add(K.layers.Dense(units=70,  activation='relu'))
                model.add(K.layers.Dense(units=50,  activation='relu'))
                model.add(K.layers.Dense(units=30, activation='relu'))
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            if m==25:
                b_size = 100
                max_epochs = 100
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=160, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dropout(0.3))
                model.add(K.layers.Dense(units=133, activation='relu'))
                model.add(K.layers.Dropout(0.24))
                model.add(K.layers.Dense(units=106, activation='relu'))
                model.add(K.layers.Dropout(0.18))
                model.add(K.layers.Dense(units=79, activation='relu'))
                model.add(K.layers.Dropout(0.12))
                model.add(K.layers.Dense(units=52, activation='relu'))
                model.add(K.layers.Dropout(0.06))
                model.add(K.layers.Dense(units=25, activation='relu'))
                model.add(K.layers.Dropout(0.02))
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            if m==30:
                b_size = 100
                max_epochs = 100
                clear_session()
                model = K.models.Sequential()
                model.add(K.layers.Dense(units=250, input_dim=train_x_zcore.shape[1],  activation='relu'))
                model.add(K.layers.Dropout(0.5))
                model.add(K.layers.Dense(units=138, activation='relu'))
                model.add(K.layers.Dropout(0.25))
                model.add(K.layers.Dense(units=25, activation='relu'))
                model.add(K.layers.Dropout(0.0))
                model.add(K.layers.Dense(units=test_y.shape[1],  activation='softmax'))
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            h = model.fit(train_x_zcore, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
            actual = np.argmax(test_y.values, axis=-1)
            predicted = np.argmax(model.predict(test_x_zcore), axis=-1)
            acc = accuracy_score(actual, predicted)
            recall = recall_score(actual, predicted, average='weighted')
            precision = precision_score(actual,predicted,average='weighted')
            F1_score=f1_score(actual,predicted,average="micro")
            #f = open("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/finall_result/"+"new_"+str(i)+"_"+str(m)+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
            f = open("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/reduplicate_7_15/"+"RF_1_"+str(m)+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
            f.write("test_acc: "+str(acc)) #将字符串写入文件中
            f.write("\n")   #换行 
            f.write("test_recall: "+str(recall)) #将字符串写入文件中
            f.write("\n")   #换行 
            f.write("test_precision: "+str(precision)) #将字符串写入文件中
            f.write("\n")   #换行 
            f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
            f.write("\n")   #换行
            other_data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/new_random_1.csv")
            other_data1 = copy.deepcopy(other_data)
            other_data1.columns = other_data1.columns+".1"
            other_data = pd.concat([other_data,other_data1],axis=1)
            other_data = other_data[train_x_zcore.columns]
            other_data_zscore = (other_data-train_x.mean(axis=0))/train_x.std(axis=0)
            predict = np.argmax(model.predict(other_data_zscore), axis=-1)
            predict2 = []
            for h in predict: 
                predict2.append(target_names[h])
            other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/new_random_1_cancertype.csv")
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