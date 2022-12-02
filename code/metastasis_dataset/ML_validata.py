import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix,precision_score, recall_score, matthews_corrcoef,cohen_kappa_score
from sklearn.metrics import roc_auc_score,roc_curve,auc


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

def cal_metrics(confusion_matrix):
    n_classes = confusion_matrix.shape[0]
    ppv1 = []
    npv1 = []
    sen1 = []
    spe1 = []
    for i in range(n_classes):
        # 逐步获取 真阳，假阳，真阴，假阴四个指标，并计算三个参数
        ALL = np.sum(confusion_matrix)
        # 对角线上是正确预测的
        TP = confusion_matrix[i, i]
        # 列加和减去正确预测是该类的假阳
        FP = np.sum(confusion_matrix[:, i]) - TP
        # 行加和减去正确预测是该类的假阴
        FN = np.sum(confusion_matrix[i, :]) - TP
        # 全部减去前面三个就是真阴
        TN = ALL - TP - FP - FN
        ppv1.append(TP/(TP+FP))
        npv1.append(TN/(TN+FN))
        sen1.append(TP/(TP+FN))
        spe1.append(TN/(TN+FP))
    return ppv1, npv1, sen1, spe1

# train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote_11_22.csv")
# train_y = np.argmax(train_y.values,axis=-1)
# y_test = np.argmax(test_y.values,axis=-1)

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_11_23.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
train_y = data['project_id'].apply(lambda x: Class_dict[x])
train_x = data.iloc[:,1:5078]

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_11_23.csv")
test_y = data['project_id'].apply(lambda x: Class_dict[x])
test_x = data.iloc[:,1:5078]

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier

RF_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/RF_feature_11_23.csv")
Extra_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Extratree_feature_11_23.csv")
xgboost_feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance/Xgboost_feature_11_23.csv")

#data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote_11_22.csv")
target_names = list(Class_dict.keys())

validata_metastasis = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata/metastatic.csv")
n = validata_metastasis.iloc[:,0]!="TCGA-ESCA"
validata_metastasis = validata_metastasis.loc[n,]
metastatic_validata_really = validata_metastasis.iloc[:,0]
#validata_Recurren = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata/recurren.csv")
#recurren_validata_really = validata_Recurren.iloc[:,0]
##################################### RF validata
geo_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/GEO_metasta_validata_cancertype.csv")
geo_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/GEO_metasta_validata.csv")
metastatic_validata_really = pd.DataFrame(metastatic_validata_really)
geo_metasta_validata_really.columns = metastatic_validata_really.columns
RF_metastatic_validata_really = pd.concat([metastatic_validata_really,geo_metasta_validata_really])
#############################    Extratree validata
EXtra_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/Extratree_GEO_metasta_validata_cancertype.csv")
Extra_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/Extratree_GEO_metasta_validata.csv")
EXtra_metasta_validata_really.columns = metastatic_validata_really.columns
EXtra_metasta_validata_really = pd.concat([metastatic_validata_really,EXtra_metasta_validata_really])
#############################    xgboost validata
xgboost_metasta_validata_really = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/xgboost_GEO_metasta_validata_cancertype.csv")
xgboost_metasta_validata = pd.read_csv("/public/slst/home/ningwei/methylation/data/GEO_validata/xgboost_GEO_metasta_validata.csv")
xgboost_metasta_validata_really.columns = metastatic_validata_really.columns
xgboost_metasta_validata_really = pd.concat([metastatic_validata_really,xgboost_metasta_validata_really])


for i in [30]:
    # 随机森林
    new_train_x = train_x[RF_feature.iloc[0:i,0]]
    new_test_x = test_x[RF_feature.iloc[0:i,0]]
    # train_x1_zcore=(new_train_x-new_train_x.mean(axis=0))/new_train_x.std(axis=0)
    # value_x_zcore=(new_test_x-new_train_x.mean(axis=0))/new_train_x.std(axis=0)
    clf2 = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
    #scores2 = cross_val_score(clf2, X, y)
    #print(scores2.mean())
    clf2.fit(new_train_x,train_y)
    ################### mestasis
    RF_data = validata_metastasis[RF_feature.iloc[0:i,0]]
    RF_data1 = geo_metasta_validata[RF_feature.iloc[0:30,0]] ### GEO-metastasic
    RF_data = pd.concat([RF_data,RF_data1])
    R_predict = clf2.predict_proba(RF_data)
    R_predict = pd.DataFrame(R_predict)
    R_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/RF_metasta.csv",index=0)
    R_predict = clf2.predict(RF_data)
    predict2 = []
    for j in R_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(RF_metastatic_validata_really,predict2)
    F1_score = f1_score(RF_metastatic_validata_really,predict2,average='weighted')
    recall = recall_score(RF_metastatic_validata_really, predict2, average='weighted')
    precision = precision_score(RF_metastatic_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(RF_metastatic_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(RF_metastatic_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"RF_metastasis_11_24"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(RF_metastatic_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/RF_metastasis_confus_11_24.csv")
    ################### recurren
    RF_data = validata_Recurren[RF_feature.iloc[0:i,0]]
    R_predict = clf2.predict(RF_data)
    predict2 = []
    for j in R_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(recurren_validata_really,predict2)
    F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
    recall = recall_score(recurren_validata_really, predict2, average='weighted')
    precision = precision_score(recurren_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(recurren_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"RF_recurren_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/RF_recurren_confus_11_23.csv")
    ############ Extratree

    new_train_x = train_x[Extra_feature.iloc[0:i,0]]
    new_test_x = test_x[Extra_feature.iloc[0:i,0]]
    clf_Extra = ExtraTreesClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
    # scores3 = cross_val_score(clf3, X, y)
    # print(scores3.mean())
    clf_Extra.fit(new_train_x,train_y)
    ################### mestasis
    Extratree_data = validata_metastasis[Extra_feature.iloc[0:i,0]]
    Extratree_data1 = Extra_metasta_validata[Extra_feature.iloc[0:30,0]] ### GEO-metastasic
    Extratree_data = pd.concat([Extratree_data,Extratree_data1])
    E_predict = clf_Extra.predict_proba(Extratree_data)
    E_predict = pd.DataFrame(E_predict)
    E_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/Extratree_metasta.csv",index=0)
    E_predict = clf_Extra.predict(Extratree_data)
    predict2 = []
    for j in E_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(EXtra_metasta_validata_really,predict2)
    F1_score = f1_score(EXtra_metasta_validata_really,predict2,average='weighted')
    recall = recall_score(EXtra_metasta_validata_really, predict2, average='weighted')
    precision = precision_score(EXtra_metasta_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(EXtra_metasta_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(EXtra_metasta_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"Extra_metastasis_11_24"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(EXtra_metasta_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/Extra_metastasis_confus_11_24.csv")
    ################### recurren
    Extratree_data = validata_Recurren[Extra_feature.iloc[0:i,0]]
    E_predict = clf_Extra.predict(Extratree_data)
    predict2 = []
    for j in E_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(recurren_validata_really,predict2)
    F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
    recall = recall_score(recurren_validata_really, predict2, average='weighted')
    precision = precision_score(recurren_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(recurren_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"Extra_recurren_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/Extra_recurren_confus_11_23.csv")
    ##################  XGBOOST
    import xgboost as xgb
    new_train_x = train_x[xgboost_feature.iloc[0:i,0]]
    new_test_x = test_x[xgboost_feature.iloc[0:i,0]]
    model = xgb.XGBClassifier(max_depth=5,learning_rate=0.1,n_estimators=160,silent=True,objective='multi:softmax')
    model.fit(new_train_x,train_y)
    ################### mestasis
    xgboost_data = validata_metastasis[xgboost_feature.iloc[0:i,0]]
    xgboost_data1 = xgboost_metasta_validata[xgboost_feature.iloc[0:30,0]] ### GEO-metastasic
    xgboost_data = pd.concat([xgboost_data,xgboost_data1])
    X_predict = model.predict_proba(xgboost_data)
    X_predict = pd.DataFrame(X_predict)
    X_predict.to_csv("/public/slst/home/ningwei/methylation/process_data/ROC_data/xgboost_metasta.csv",index=0)
    X_predict = model.predict(xgboost_data)
    predict2 = []
    for j in X_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(xgboost_metasta_validata_really,predict2)
    F1_score = f1_score(xgboost_metasta_validata_really,predict2,average='weighted')
    recall = recall_score(xgboost_metasta_validata_really, predict2, average='weighted')
    precision = precision_score(xgboost_metasta_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(xgboost_metasta_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(xgboost_metasta_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"xgboost_metastasis_11_24"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(xgboost_metasta_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/xgboost_metastasis_confus_11_24.csv")
    ################### recurren
    xgboost_data = validata_Recurren[xgboost_feature.iloc[0:i,0]]
    X_predict = model.predict(xgboost_data)
    predict2 = []
    for j in X_predict:
        predict2.append(target_names[j])
    accuracy = accuracy_score(recurren_validata_really,predict2)
    F1_score = f1_score(recurren_validata_really,predict2,average='weighted')
    recall = recall_score(recurren_validata_really, predict2, average='weighted')
    precision = precision_score(recurren_validata_really, predict2,average='weighted')
    matthews = matthews_corrcoef(recurren_validata_really, predict2,  sample_weight=None)
    kappa = cohen_kappa_score(recurren_validata_really, predict2)
    f = open("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/"+"xgboost_recurren_11_23"+".txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
    f.write("validata_acc: "+str(accuracy)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_recall: "+str(recall)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_precision: "+str(precision)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_F1: "+str(F1_score)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_matthew: "+str(matthews)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_kappa: "+str(kappa)) #将字符串写入文件中
    f.write("\n")   #换行
    predict2 = pd.DataFrame(predict2)
    confusion_mat = confusion_matrix(recurren_validata_really,predict2,labels=target_names)
    ppv, npv, sen, spe = cal_metrics(confusion_mat)
    f.write("validata_ppv: "+str(ppv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_npv: "+str(npv)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_sen: "+str(sen)) #将字符串写入文件中
    f.write("\n")   #换行
    f.write("validata_spe: "+str(spe)) #将字符串写入文件中
    f.write("\n")   #换行
    f.close()
    confusion_mat1 = pd.DataFrame(confusion_mat)
    confusion_mat1.to_csv("/public/slst/home/ningwei/methylation/process_data/TCGA_validata_result/xgboost_recurren_confus_11_23.csv")
