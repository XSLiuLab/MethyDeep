import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from collections import Counter
import sys
# def load_data(CSV_FILE_PATH):
#     IRIS = pd.read_csv(CSV_FILE_PATH)
#     target_var = 'project_id'  # 目标变量
#     # 数据集的特征
#     features = list(IRIS.columns)
#     features.remove(target_var)
#     # 目标变量的类别
#     Class = IRIS[target_var].unique()
#     # 目标变量的类别字典
#     Class_dict = dict(zip(Class, range(len(Class))))
#     # 增加一列target, 将目标变量进行编码
#     IRIS['target'] = IRIS[target_var].apply(lambda x: Class_dict[x])
#     # 对目标变量进行0-1编码(One-hot Encoding)
#     lb = LabelBinarizer()
#     lb.fit(list(Class_dict.values()))
#     transformed_labels = lb.transform(IRIS['target'])
#     y_bin_labels = []  # 对多分类进行0-1编码的变量
#     for i in range(transformed_labels.shape[1]):
#         y_bin_labels.append('y' + str(i))
#         IRIS['y' + str(i)] = transformed_labels[:, i]
#     # 将数据集分为训练集和测试集
#     train_x, test_x, train_y, test_y = train_test_split(IRIS[features], IRIS[y_bin_labels], stratify=IRIS['project_id'],\
#                                                         train_size=0.8, test_size=0.2, random_state=1)
#     return train_x, test_x, train_y, test_y

# train_x, test_x, train_y, test_y = load_data("/public/slst/home/ningwei/methylation/process_data/train_data/data_smote_11_22.csv")
# train_y = np.argmax(train_y.values,axis=-1)
# test_y = np.argmax(test_y.values,axis=-1)

m = sys.argv[1]
print(m)
m=int(m)


data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv")
Class = data['project_id'].unique()
Class_dict = dict(zip(Class, range(len(Class))))
train_y = data['project_id'].apply(lambda x: Class_dict[x])
train_x = data.iloc[:,1:data.shape[1]]

data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_6.csv")
test_y = data['project_id'].apply(lambda x: Class_dict[x])
test_x = data.iloc[:,1:data.shape[1]]

feature = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/cg_frequency.csv")
new_feature = feature.loc[feature['Freq']>=m,]

train_x = train_x[new_feature['all_cg']]
test_x = test_x[new_feature['all_cg']]

###################################
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier

## 随机森林
clf2 = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
#scores2 = cross_val_score(clf2, X, y)
#print(scores2.mean())
clf2.fit(train_x,train_y)
random_forest_importance = clf2.feature_importances_
random_forest_feature_importance=[(feature,round(importance,8)) 
                                  for feature, importance in zip(train_x.columns,random_forest_importance)]
random_forest_feature_importance=sorted(random_forest_feature_importance,key=lambda x:x[1],reverse=True)
random_forest_feature_importance = pd.DataFrame(random_forest_feature_importance)
random_forest_feature_importance.to_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/RF_feature"+str(m)+".csv",index=0)
############# Extratree
clf_Extra = ExtraTreesClassifier(n_estimators=10, max_depth=None,min_samples_split=2, random_state=0)
# scores3 = cross_val_score(clf3, X, y)
# print(scores3.mean())
clf_Extra.fit(train_x,train_y)
Extratree_importance = clf_Extra.feature_importances_
Extratree_feature_importance=[(feature,round(importance,8)) 
                                  for feature, importance in zip(train_x.columns,Extratree_importance)]
Extratree_feature_importance=sorted(Extratree_feature_importance,key=lambda x:x[1],reverse=True)
Extratree_feature_importance = pd.DataFrame(Extratree_feature_importance)
Extratree_feature_importance.to_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/Extratree_feature"+str(m)+".csv",index=0)
##
###################  XGBOOST
import xgboost as xgb
model = xgb.XGBClassifier(max_depth=5,learning_rate=0.1,n_estimators=160,silent=True,objective='multi:softmax')
model.fit(train_x,train_y)
xgboost_importance = model.feature_importances_
xgboost_feature_importance=[(feature,round(importance,8)) 
                                  for feature, importance in zip(train_x.columns,xgboost_importance)]
xgboost_feature_importance=sorted(xgboost_feature_importance,key=lambda x:x[1],reverse=True)
xgboost_feature_importance = pd.DataFrame(xgboost_feature_importance)
xgboost_feature_importance.to_csv("/public/slst/home/ningwei/methylation/process_data/feature_importance_12_6/Xgboost_feature"+str(m)+".csv",index=0)

