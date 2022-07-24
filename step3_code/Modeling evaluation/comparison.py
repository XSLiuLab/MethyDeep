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

def load_data(CSV_FILE_PATH):
    IRIS = pd.read_csv(CSV_FILE_PATH)
    target_var = 'CancerTypes'  # 目标变量
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
    train_x, test_x, train_y, test_y = train_test_split(IRIS[features], IRIS[y_bin_labels], stratify=IRIS['CancerTypes'],\
                                                        train_size=0.8, test_size=0.2, random_state=1)
    return train_x, test_x, train_y, test_y

m=["EXtratree","lasso","RF"]

for i in range(0,3) :
    for j in [5,10,15,20,25,30]:
        # if i ==0 :
        #     train_x, test_x, train_y, test_y = load_data("/slst/home/ningwei/TCGA/Methylation450K/9_17_50cg_methydata_cancertype.csv")
        # else:
        #train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/TCGA/Methylation450K/Three_method/data/"+m[i]+".csv")
        train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/TCGA/Methylation450K/data/champ_data/"+m[i]+"_50.csv")
        train_x = train_x.iloc[:,0:j]
        test_x = test_x.iloc[:,0:j]
        print(j)
        print(train_x.shape)
        #train_x1, value_x, train_y1, value_y = train_test_split(train_x, train_y,train_size=0.8, test_size=0.2, random_state=1)
        ###############################################################  对数据进行标准化
        train_x_zcore=(train_x-train_x.mean(axis=0))/train_x.std(axis=0)
        test_x_zcore=(test_x-train_x.mean(axis=0))/train_x.std(axis=0)
        ############################## EXtratree
        if i == 0 and j==5 :
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=500, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=44, activation='relu'))      
            model.add(K.layers.Dense(units=38, activation='relu'))  
            model.add(K.layers.Dense(units=32, activation='relu')) 
            model.add(K.layers.Dense(units=25, activation='relu'))  
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 0 and j==10 :
            b_size = 100
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=100, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=63, activation='relu'))      
            model.add(K.layers.Dense(units=25, activation='relu'))  
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 0 and j==15 :
            b_size = 100
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=130, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=80, activation='relu'))      
            model.add(K.layers.Dense(units=30, activation='relu'))  
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 0 and j==20 :
            b_size = 50
            max_epochs = 50
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=150, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=120, activation='relu'))   
            model.add(K.layers.Dropout(0.23))
            model.add(K.layers.Dense(units=90, activation='relu'))      
            model.add(K.layers.Dropout(0.16))
            model.add(K.layers.Dense(units=60, activation='relu'))   
            model.add(K.layers.Dropout(0.09))
            model.add(K.layers.Dense(units=30, activation='relu'))  
            model.add(K.layers.Dropout(0.02))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 0 and j==25 :
            b_size = 100
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=200, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=144, activation='relu'))   
            model.add(K.layers.Dropout(0.2))
            model.add(K.layers.Dense(units=87, activation='relu'))      
            model.add(K.layers.Dropout(0.1))
            model.add(K.layers.Dense(units=30, activation='relu'))  
            model.add(K.layers.Dropout(0.0))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 0 and j==30 :
            b_size = 100
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=230, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.5))
            model.add(K.layers.Dense(units=179, activation='relu'))   
            model.add(K.layers.Dropout(0.38))
            model.add(K.layers.Dense(units=128, activation='relu'))      
            model.add(K.layers.Dropout(0.26))
            model.add(K.layers.Dense(units=77, activation='relu'))      
            model.add(K.layers.Dropout(0.14))
            model.add(K.layers.Dense(units=25, activation='relu'))  
            model.add(K.layers.Dropout(0.02))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    
        ####################### LASSO
        if i ==1 and j==5 :
            b_size = 80
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=40, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=33, activation='relu'))
            model.add(K.layers.Dense(units=25, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i ==1 and j==10 :
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=60, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=43, activation='relu'))
            model.add(K.layers.Dense(units=25, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i ==1 and j==15 :
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=150, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=88, activation='relu'))
            model.add(K.layers.Dropout(0.15))
            model.add(K.layers.Dense(units=25, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i ==1 and j==20 :
            b_size = 80
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=200, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=144, activation='relu'))
            model.add(K.layers.Dropout(0.2))
            model.add(K.layers.Dense(units=87, activation='relu'))
            model.add(K.layers.Dropout(0.1))
            model.add(K.layers.Dense(units=30, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 1 and j==25 :
            b_size = 80
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=250, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.5))
            model.add(K.layers.Dense(units=140, activation='relu'))   
            model.add(K.layers.Dropout(0.25))
            model.add(K.layers.Dense(units=30, activation='relu'))  
            model.add(K.layers.Dropout(0.0))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 1 and j==30 :
            b_size = 100
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=300, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.5))
            model.add(K.layers.Dense(units=163, activation='relu'))   
            model.add(K.layers.Dropout(0.25))
            model.add(K.layers.Dense(units=25, activation='relu'))  
            model.add(K.layers.Dropout(0.0))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        ######################## RF
        if i == 2 and j==5:
            b_size = 80
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=40, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dense(units=37, activation='relu'))
            model.add(K.layers.Dense(units=34, activation='relu'))
            model.add(K.layers.Dense(units=30, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 2 and j==10:
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=100, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=77, activation='relu'))
            model.add(K.layers.Dropout(0.2))
            model.add(K.layers.Dense(units=54, activation='relu'))
            model.add(K.layers.Dropout(0.1))
            model.add(K.layers.Dense(units=30, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 2 and j==15:
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=110, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=84, activation='relu'))
            model.add(K.layers.Dropout(0.2))
            model.add(K.layers.Dense(units=57, activation='relu'))
            model.add(K.layers.Dropout(0.1))
            model.add(K.layers.Dense(units=30, activation='relu'))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 2 and j==20:
            b_size = 80
            max_epochs = 50
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=130, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.3))
            model.add(K.layers.Dense(units=104, activation='relu'))
            model.add(K.layers.Dropout(0.23))
            model.add(K.layers.Dense(units=78, activation='relu'))
            model.add(K.layers.Dropout(0.16))
            model.add(K.layers.Dense(units=52, activation='relu'))
            model.add(K.layers.Dropout(0.09))
            model.add(K.layers.Dense(units=25, activation='relu'))
            model.add(K.layers.Dropout(0.02))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 2 and j==25:
            b_size = 50
            max_epochs = 100
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=200, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.5))
            model.add(K.layers.Dense(units=142, activation='relu'))
            model.add(K.layers.Dropout(0.33))
            model.add(K.layers.Dense(units=84, activation='relu'))
            model.add(K.layers.Dropout(0.16))
            model.add(K.layers.Dense(units=25, activation='relu'))
            model.add(K.layers.Dropout(0.0))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        if i == 2 and j==30:
            b_size = 100
            max_epochs = 50
            clear_session()
            model = K.models.Sequential()
            model.add(K.layers.Dense(units=280, input_dim=train_x_zcore.shape[1],  activation='relu'))
            model.add(K.layers.Dropout(0.5))
            model.add(K.layers.Dense(units=197, activation='relu'))
            model.add(K.layers.Dropout(0.33))
            model.add(K.layers.Dense(units=114, activation='relu'))
            model.add(K.layers.Dropout(0.16))
            model.add(K.layers.Dense(units=31, activation='relu'))
            model.add(K.layers.Dropout(0.0))
            model.add(K.layers.Dense(units=24,  activation='softmax'))
            model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        h = model.fit(train_x_zcore, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
        actual = np.argmax(test_y.values, axis=-1)
        predicted = np.argmax(model.predict(test_x_zcore), axis=-1)
        acc = accuracy_score(actual, predicted)
        recall = recall_score(actual, predicted, average='weighted')
        precision = precision_score(actual,predicted,average='weighted')
        F1_score = f1_score(actual,predicted,average="micro")
        f = open("/public/slst/home/ningwei/TCGA/Methylation450K/data/result/multiple_"+m[i]+"_"+ str(j) +".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
        f.write("test_acc: "+str(acc)) #将字符串写入文件中
        f.write("\n")   #换行 
        f.write("test_recall: "+str(recall)) #将字符串写入文件中
        f.write("\n")   #换行 
        f.write("test_precision: "+str(precision)) #将字符串写入文件中
        f.write("\n")   #换行
        f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
        f.write("\n")   #换行
        #data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/Three_method/data/"+m[i]+".csv")
        data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/data/champ_data/"+m[i]+"_50.csv")
        target_names = data.iloc[:,0].unique()
        other_data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/data/validata/"+m[i]+".csv")
        other_data = other_data[train_x_zcore.columns]
        other_data_zscore = (other_data-train_x.mean(axis=0))/train_x.std(axis=0)
        predict = np.argmax(model.predict(other_data_zscore), axis=-1)
        predict2 = []
        for j in predict: 
            predict2.append(target_names[j])
        #other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[i]+"_cancertype.csv")
        other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/data/validata/"+m[i]+"_cancertype.csv")
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