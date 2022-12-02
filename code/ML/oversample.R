######################  Oversampling and splitting data
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import RandomOverSampler
#Define SMOTE model, random_ State is equivalent to random number seed
data = pd.read_csv("/public/slst/home/ningwei/methylation/process_data/train_data/data_11_19.csv")
y = data.iloc[:,0]
x = data.iloc[:,1:5078]
print(Counter(y))

y = y.replace('TCGA-LGG', 'TCGA-GBM') ### merge LGG_GBM
y = y.replace('TCGA-READ', 'TCGA-COAD') ### merge READ COAD

x = x.loc[y!="TCGA-STAD",]  #### delete STAD
x = x.loc[y!="TCGA-ESCA",]  #### delete ESCA

y = y.loc[y!="TCGA-STAD",]  #### delete STAD
y = y.loc[y!="TCGA-ESCA",]  #### delete ESCA

train_x, test_x, train_y, test_y = train_test_split(x, y, stratify=y,\
                                                        train_size=0.8, test_size=0.2, random_state=1)

# smo = SMOTE(random_state=42)
# X_smo, y_smo = smo.fit_sample(train_x, train_y)
# print(Counter(y_smo))

ros = RandomOverSampler(random_state=0)
X_resampled, y_resampled = ros.fit_resample(train_x, train_y)
counter_resampled = Counter(y_resampled)
print(counter_resampled)

data = pd.concat([y_resampled,X_resampled],axis=1)
data.to_csv("/public/slst/home/ningwei/methylation/process_data/train_data/train_data_smote_11_23.csv",index=0)

data = pd.concat([test_y,test_x],axis=1)
data.to_csv("/public/slst/home/ningwei/methylation/process_data/train_data/test_data_smote_11_23.csv",index=0)