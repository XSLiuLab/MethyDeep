######################  Oversampling and splitting data
from imblearn.over_sampling import RandomOverSampler
#Define SMOTE model, random_ State is equivalent to random number seed
data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/data_12_5.csv")
y = data.iloc[:,0]
x = data.iloc[:,1:data.shape[1]]
print(Counter(y))

y = y.replace('TCGA-LGG', 'TCGA-GBM') ### 合并LGG_GBM
y = y.replace('TCGA-READ', 'TCGA-COAD') ### 合并 READ COAD

x = x.loc[y!="TCGA-STAD",]  #### 删除 STAD
y = y.loc[y!="TCGA-STAD",]  #### 删除 STAD

x = x.loc[y!="TCGA-ESCA",]  #### 删除 ESCA
y = y.loc[y!="TCGA-ESCA",]  #### 删除 ESCA

train_x, test_x, train_y, test_y = train_test_split(x, y, stratify=y,\
                                                        train_size=0.8, test_size=0.2, random_state=1)



ros = RandomOverSampler(random_state=0)
X_resampled, y_resampled = ros.fit_resample(train_x, train_y)
counter_resampled = Counter(y_resampled)
print(counter_resampled)

data = pd.concat([y_resampled,X_resampled],axis=1)
data.to_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_12_6.csv",index=0)

data = pd.concat([test_y,test_x],axis=1)
data.to_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_12_6.csv",index=0)