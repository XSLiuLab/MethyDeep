import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score


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



from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier

RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/RF_feature_11_23.csv")
Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Extratree_feature_11_23.csv")
xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")


for i in [5,10,15,20,25,30,35,40,45,50]:
    ## RF
    new_train_x = train_x[RF_feature.iloc[0:i,0]]
    new_test_x = test_x[RF_feature.iloc[0:i,0]]
    clf2 = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
    #scores2 = cross_val_score(clf2, X, y)
    #print(scores2.mean())
    clf2.fit(new_train_x,train_y)
    R_predict = clf2.predict(new_test_x)
    accuracy = accuracy_score(test_y,R_predict)
    F1_score = f1_score(test_y,R_predict,average='macro')
    recall = recall_score(test_y, R_predict, average='weighted')
    precision = precision_score(test_y, R_predict,average='weighted')
    f = open("/public/slst/home/ningwei/methylation/process_data/ML_different_number/"+"随机森林"+ str(i)+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(accuracy)) 
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) 
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) 
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) 
    f.write("\n")   #换行
    f.close()
    
    ############# Extratree
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
    f = open("/public/slst/home/ningwei/methylation/process_data/ML_different_number/"+"ExtraTree"+str(i)+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(accuracy)) 
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) 
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) 
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) 
    f.write("\n")   #换行
    f.close()
    ###################  XGBOOST
    import xgboost as xgb
    new_train_x = train_x[xgboost_feature.iloc[0:i,0]]
    new_test_x = test_x[xgboost_feature.iloc[0:i,0]]
    model = xgb.XGBClassifier(max_depth=5,learning_rate=0.1,n_estimators=160,silent=True,objective='multi:softmax')
    model.fit(new_train_x,train_y)
    X_pred = model.predict(new_test_x)
    accuracy = accuracy_score(test_y,X_pred)
    F1_score = f1_score(test_y, X_pred,average='macro')
    recall = recall_score(test_y, X_pred, average='weighted')
    precision = precision_score(test_y, X_pred,average='weighted')
    f = open("/public/slst/home/ningwei/methylation/process_data/ML_different_number/"+"XGBOOST"+str(i)+"_11_23.txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("test_acc: "+str(accuracy)) 
    f.write("\n")   #换行 
    f.write("test_recall: "+str(recall)) 
    f.write("\n")   #换行 
    f.write("test_precision: "+str(precision)) 
    f.write("\n")   #换行
    f.write("test_F1: "+str(F1_score)) 
    f.write("\n")   #换行
    f.close()