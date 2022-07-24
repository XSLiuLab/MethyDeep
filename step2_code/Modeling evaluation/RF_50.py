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

for j in range(0,100):
    for l in [0,3,8,13,18,23,28]:
        for m in [5,10,15,20,25,30]:
            train_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_x_"+str(l)+"_"+str(m)+".csv")
            train_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_y_"+str(l)+"_"+str(m)+".csv")
            test_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_test_x_"+str(l)+"_"+str(m)+".csv")
            test_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_test_y_"+str(l)+"_"+str(m)+".csv")
            data = pd.read_csv("~/TCGA/Methylation450K/data_6_15/"+gred_list[l])
            target_names = data.iloc[:,0].unique()
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
            # if m==35:
            #     b_size = 80
            #     max_epochs = 100
            #     clear_session()
            #     model = K.models.Sequential()
            #     model.add(K.layers.Dense(units=280, input_dim=train_x_zcore.shape[1],  activation='relu'))
            #     model.add(K.layers.Dropout(0.3))
            #     model.add(K.layers.Dense(units=153, activation='relu'))
            #     model.add(K.layers.Dropout(0.15))
            #     model.add(K.layers.Dense(units=25, activation='relu'))
            #     model.add(K.layers.Dropout(0.0))
            #     model.add(K.layers.Dense(units=24,  activation='softmax'))
            #     model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            # if m==40:
            #     b_size = 80
            #     max_epochs = 50
            #     clear_session()
            #     model = K.models.Sequential()
            #     model.add(K.layers.Dense(units=250, input_dim=train_x_zcore.shape[1],  activation='relu'))
            #     model.add(K.layers.Dropout(0.3))
            #     model.add(K.layers.Dense(units=194, activation='relu'))
            #     model.add(K.layers.Dropout(0.23))
            #     model.add(K.layers.Dense(units=138, activation='relu'))
            #     model.add(K.layers.Dropout(0.16))
            #     model.add(K.layers.Dense(units=82, activation='relu'))
            #     model.add(K.layers.Dropout(0.09))
            #     model.add(K.layers.Dense(units=25, activation='relu'))
            #     model.add(K.layers.Dropout(0.02))
            #     model.add(K.layers.Dense(units=24,  activation='softmax'))
            #     model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            # if m==45:
            #     b_size = 100
            #     max_epochs = 100
            #     clear_session()
            #     model = K.models.Sequential()
            #     model.add(K.layers.Dense(units=310, input_dim=train_x_zcore.shape[1],  activation='relu'))
            #     model.add(K.layers.Dropout(0.5))
            #     model.add(K.layers.Dense(units=215, activation='relu'))
            #     model.add(K.layers.Dropout(0.33))
            #     model.add(K.layers.Dense(units=120, activation='relu'))
            #     model.add(K.layers.Dropout(0.16))
            #     model.add(K.layers.Dense(units=25, activation='relu'))
            #     model.add(K.layers.Dropout(0.0))
            #     model.add(K.layers.Dense(units=24,  activation='softmax'))
            #     model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            # if m==50:
            #     b_size = 100
            #     max_epochs = 50
            #     clear_session()
            #     model = K.models.Sequential()
            #     model.add(K.layers.Dense(units=280, input_dim=train_x_zcore.shape[1],  activation='relu'))
            #     model.add(K.layers.Dropout(0.3))
            #     model.add(K.layers.Dense(units=153, activation='relu'))
            #     model.add(K.layers.Dropout(0.15))
            #     model.add(K.layers.Dense(units=25, activation='relu'))
            #     model.add(K.layers.Dense(units=24,  activation='softmax'))
            #     model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
            h = model.fit(train_x_zcore, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
            actual = np.argmax(test_y.values, axis=-1)
            predicted = np.argmax(model.predict(test_x_zcore), axis=-1)
            acc = accuracy_score(actual, predicted)
            recall = recall_score(actual, predicted, average='weighted')
            precision = precision_score(actual,predicted,average='weighted')
            if l==0:
                q=2
            else:
                q=l+2
            f = open("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/reduplicate_7_15/"+"RF_"+str(q)+"_"+str(m)+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
            f.write("test_acc: "+str(acc)) #将字符串写入文件中
            f.write("\n")   #换行 
            f.write("test_recall: "+str(recall)) #将字符串写入文件中
            f.write("\n")   #换行 
            f.write("test_precision: "+str(precision)) #将字符串写入文件中
            f.write("\n")   #换行 
            if q != 2:
                q=q-2
            other_data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/random_"+str(q)+".csv")
            other_data1 = copy.deepcopy(other_data)
            other_data1.columns = other_data1.columns+".1"
            other_data2 = copy.deepcopy(other_data)
            other_data2.columns = other_data2.columns+".2"
            other_data3 = copy.deepcopy(other_data)
            other_data3.columns = other_data3.columns+".3"
            other_data4 = copy.deepcopy(other_data)
            other_data4.columns = other_data4.columns+".4"
            other_data5 = copy.deepcopy(other_data)
            other_data5.columns = other_data5.columns+".5"
            other_data6 = copy.deepcopy(other_data)
            other_data6.columns = other_data6.columns+".6"
            other_data7 = copy.deepcopy(other_data)
            other_data7.columns = other_data7.columns+".7"
            other_data8 = copy.deepcopy(other_data)
            other_data8.columns = other_data8.columns+".8"
            other_data9 = copy.deepcopy(other_data)
            other_data9.columns = other_data9.columns+".9"
            other_data = pd.concat([other_data,other_data1,other_data2,other_data3,other_data4,other_data5,other_data6,other_data7,other_data8,other_data9],axis=1)
            other_data = other_data[train_x_zcore.columns]
            other_data_zscore = (other_data-train_x.mean(axis=0))/train_x.std(axis=0)
            predict = np.argmax(model.predict(other_data_zscore), axis=-1)
            predict2 = []
            for h in predict: 
                predict2.append(target_names[h])
            other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/validata_mutlita/random_"+str(q)+"_cancertype.csv")
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
