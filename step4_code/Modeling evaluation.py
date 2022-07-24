import numpy as np
from sklearn.linear_model import LogisticRegression as LR
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
import pandas as pd
from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import train_test_split
from xgboost.sklearn import XGBClassifier
import sys
from tensorflow import keras as K
import tensorflow as tf
from tensorflow.keras import regularizers
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedShuffleSplit
import matplotlib.pyplot as plt
from tensorflow.keras.backend import clear_session


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
m=["Extratree","lasso","RF"]
for h in range(0,3) :
    train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/TCGA/Methylation450K/Three_method/data/"+m[h]+".csv") #
    train_x = train_x.iloc[:,0:30]
    test_x = test_x.iloc[:,0:30]
    train_x_zcore=(train_x-train_x.mean(axis=0))/train_x.std(axis=0)
    test_x_zcore=(test_x-train_x.mean(axis=0))/train_x.std(axis=0)
    train_really = np.argmax(train_y.values, axis=-1)
    test_really = np.argmax(test_y.values, axis=-1)
    ########################## LR
    lr = LR(multi_class='multinomial',solver='lbfgs',class_weight='balanced',max_iter=1000)
    #For multiclass problems, only 'newton-cg', 'sag', 'saga' and 'lbfgs' handle multinomial loss
    lr.fit(train_x_zcore,train_really) ##拟合模型
    score = lr.score(train_x_zcore,train_really)
    print(score)
    test_score = lr.score(test_x_zcore,test_really)
    print(test_score)#0.6
    acc = accuracy_score(test_really, lr.predict(test_x_zcore))
    recall = recall_score(test_really, lr.predict(test_x_zcore), average='weighted')
    precision = precision_score(test_really, lr.predict(test_x_zcore),average='weighted')
    F1_score = f1_score(test_really, lr.predict(test_x_zcore),average="micro")
    f = open("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/result/LR__"+m[h]+"_30.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/Three_method/data/"+m[h]+".csv")
    target_names = data.iloc[:,0].unique()
    other_data = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[h]+".csv")
    other_data = other_data[train_x_zcore.columns]
    other_data_zscore = (other_data-train_x.mean(axis=0))/train_x.std(axis=0)
    predict = lr.predict(other_data_zscore)
    predict2 = []
    for i in predict: 
        predict2.append(target_names[i])
    other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[h]+"_cancertype.csv")
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
    ######################### RF
    random_forest_seed=1234
    random_forest_model=RandomForestClassifier(n_estimators=200,random_state=random_forest_seed)
    random_forest_model.fit(train_x_zcore,train_really)
    rf_predictions = random_forest_model.predict(test_x_zcore)
    acc = accuracy_score(test_really, rf_predictions)
    recall = recall_score(test_really, rf_predictions, average='weighted')
    precision = precision_score(test_really, rf_predictions,average='weighted')
    F1_score = f1_score(test_really, rf_predictions,average="micro")
    f = open("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/result/RF__"+m[h]+"_30.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    predict = random_forest_model.predict(other_data_zscore)
    predict2 = []
    for i in predict: 
        predict2.append(target_names[i])
    other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[h]+"_cancertype.csv")
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
    ####################### XGBOOTS
    xgboost_clf = XGBClassifier(min_child_weight=6,max_depth=15,objective='multi:softmax',num_class=24)
    xgboost_clf.fit(train_x_zcore, train_really)
    pre_y_test = xgboost_clf.predict(test_x_zcore)
    acc = accuracy_score(test_really, pre_y_test)
    recall = recall_score(test_really, pre_y_test, average='weighted')
    precision = precision_score(test_really, pre_y_test,average='weighted')
    F1_score = f1_score(test_really, pre_y_test,average="micro")
    f = open("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/result/Xgboots__"+m[h]+"_30.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    predict = xgboost_clf.predict(other_data_zscore)
    predict2 = []
    for i in predict: 
        predict2.append(target_names[i])
    other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[h]+"_cancertype.csv")
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
    ####################### DNN
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
    h1 = model.fit(train_x_zcore, train_y, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)##调参，节点数，层数
    actual = np.argmax(test_y.values, axis=-1)
    predicted = np.argmax(model.predict(test_x_zcore), axis=-1)
    acc = accuracy_score(actual, predicted)
    recall = recall_score(actual, predicted, average='weighted')
    precision = precision_score(actual,predicted,average='weighted')
    F1_score = f1_score(actual,predicted,average="micro")
    f = open("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/result/DNN__"+m[h]+"_30.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(acc)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    predict = np.argmax(model.predict(other_data_zscore), axis=-1)
    predict2 = []
    for j in predict: 
        predict2.append(target_names[j])
    other_really = pd.read_csv("/public/slst/home/ningwei/TCGA/Methylation450K/model_comparison/validata/"+m[h]+"_cancertype.csv")
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