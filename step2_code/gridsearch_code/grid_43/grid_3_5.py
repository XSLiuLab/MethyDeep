#########################################################  加载模块
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
from sklearn.model_selection import PredefinedSplit
import math
import pandas as pd
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
########################################################################### 构建划分函数

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

def FindLayerNodesLinear(n_layers, first_layer_nodes, last_layer_nodes):
    layers = []
    nodes_increment = (last_layer_nodes - first_layer_nodes)/ (n_layers-1)
    nodes = first_layer_nodes
    for i in range(1, n_layers+1):
        layers.append(math.ceil(nodes))
        nodes = nodes + nodes_increment
    return layers

def FinddropoutLinear(n_layers, dropout):
    layers = []
    nodes_increment = round(dropout/(n_layers-1),2)
    nodes = dropout
    for i in range(1, n_layers+1):
        layers.append(nodes)
        nodes = round(nodes - nodes_increment,2)
        if(nodes <= 0):
            nodes = 0
    
    return layers


def createmodel(n_layers, first_layer_nodes, last_layer_nodes, activation_func, loss_func,dropout):
    model = Sequential()
    n_nodes = FindLayerNodesLinear(n_layers, first_layer_nodes, last_layer_nodes)
    n_dropout = FinddropoutLinear(n_layers,dropout)
    for i in range(1, n_layers):
        if i==1:
            model.add(Dense(first_layer_nodes, input_dim=train_x.shape[1], activation=activation_func))
            model.add(K.layers.Dropout(rate=n_dropout[0]))
        else:
            model.add(Dense(n_nodes[i-1], activation=activation_func))
            model.add(K.layers.Dropout(rate=n_dropout[i-1]))
    model.add(Dense(train_y.shape[1], activation='softmax'))
    model.compile(optimizer='adam', loss=loss_func, metrics = ["accuracy"]) #note: metrics could also be 'mse'
    
    return model
##################################################################  对数据进行处理
# train_x, test_x, train_y, test_y, Class_dict = load_data("methydata_numeric.csv")
# data = pd.read_csv("/methydata_numeric.csv")
# traindata = data.iloc[train_y.index,]
# traindata.to_csv("traindata.csv",index=0)
# train_x1, value_x, train_y1, value_y, Class_dict = load_data("traindata.csv")
# train_x, test_x, train_y, test_y, Class_dict = load_data("~/TCGA/Methylation450K/12_18_5551_methydata.csv")
# data = pd.read_csv("~/TCGA/Methylation450K/12_18_5551_methydata.csv")

# train_x, test_x, train_y, test_y, Class_dict = load_data("~/TCGA/Methylation450K/886_methydata.csv")   ##### 去除高相关数据集
# data = pd.read_csv("~/TCGA/Methylation450K/886_methydata.csv")

# train_x, test_x, train_y, test_y, Class_dict = load_data("~/TCGA/Methylation450K/five_cancer/151.methydata.csv")   ##### 去除高相关数据集
# data = pd.read_csv("~/TCGA/Methylation450K/five_cancer/151.methydata.csv")
# traindata = data.iloc[train_y.index,]
# traindata.to_csv("/slst/home/ningwei/TCGA/Methylation450K/110traindata.csv",index=0)
# train_x1, value_x, train_y1, value_y, Class_dict = load_data("/slst/home/ningwei/TCGA/Methylation450K/110traindata.csv")
# train_x, test_x, train_y, test_y,Class_dict = load_data("~/TCGA/Methylation450K/five_cancer/new_data.csv") ### five cancer
# train_x1, value_x, train_y1, value_y = train_test_split(train_x, train_y,train_size=0.8, test_size=0.2, random_state=1)

#train_x, test_x, train_y, test_y,Class_dict = load_data("~/TCGA/Methylation450K/five_cancer/20_new_data.csv") ### five cancer

# train_x = pd.read_csv("~/TCGA/Methylation450K/five_cancer/train_data/20_train_x.csv")
# train_y = pd.read_csv("~/TCGA/Methylation450K/five_cancer/train_data/20_train_y.csv")
for i in [43]:
    for m in [5]:
        train_x = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_x_"+str(i)+"_"+str(m)+".csv")
        train_y = pd.read_csv("~/TCGA/Methylation450K/feature_importance_6_15/RF_gredsearch/RF_train_y_"+str(i)+"_"+str(m)+".csv")
        train_x1, value_x, train_y1, value_y = train_test_split(train_x, train_y,train_size=0.8, test_size=0.2, random_state=1)
        ###############################################################  对数据进行标准化
        train_x1_zcore=(train_x1-train_x.mean(axis=0))/train_x.std(axis=0)
        value_x_zcore=(value_x-train_x.mean(axis=0))/train_x.std(axis=0)
        ############################################################################ 构建函数
        train_val_features = np.concatenate((train_x1_zcore,value_x_zcore),axis=0)
        train_val_labels = np.concatenate((train_y1,value_y ),axis=0)
        test_fold = np.zeros(train_val_features.shape[0]) 
        test_fold[:train_x1.shape[0]] = -1  
        ps = PredefinedSplit(test_fold=test_fold)
        #####################################################################设置参数范围
        model =  KerasClassifier(build_fn=createmodel, verbose = False)  

        activation_funcs = ['sigmoid', 'relu'] 
        #activation_funcs = ['relu'] 
        loss_funcs = ['categorical_crossentropy']
        if m ==5 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [50,40,30], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m ==10 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [100,80,60,50], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 15 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [150,130,110,100,80], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 20 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [200,170,150,130,110,100], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 25 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [250,230,200,180,160,130,110,100], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 30 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [300,280,250,230,200,180,150,130,100], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 35 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [350,330,310,280,250,230,200,180,150,130], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 40 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [400,380,360,340,310,280,250,230,200,180,150], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 45 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [450,430,410,400,380,360,340,310,280,250,230,200], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        if m == 50 :
            #param_grid = dict(n_layers=[2,3,4,5,6,7], first_layer_nodes = [10000,8000,6000,5000,4000], last_layer_nodes = [60,90,120], dropout=[0,0.1,0.3,0.5,0.7,0.9], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,50], epochs = [50])
            param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [500,470,450,430,410,400,380,360,340,310,280,250,230,200], last_layer_nodes = [30,25], dropout=[0,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])
        grid = GridSearchCV(estimator = model, param_grid = param_grid,cv=ps,n_jobs=1)
        #grid = RandomizedSearchCV (estimator = model, param_grid = param_grid,cv=3,n_jobs=5)
        ################################################################ 进行训练
        grid.fit(train_val_features, train_val_labels)
        ############################################################### 输出最高精确度及相应参数
        print(grid.best_score_)
        print(grid.best_params_)
        f = open("/public/slst/home/ningwei/TCGA/Methylation450K/feature_importance_6_15/best_param/"+"RF_"+str(i)+"_"+str(m)+"txt",'a') #若文件不存在，系统自动创建。'a'表示可连续写入到文件，保留原内容，在原内容之后写入。可修改该模式（'w+','w','wb'等）
        f.write(grid.best_score_) #将字符串写入文件中
        f.write("\n")   #换行 
        f.write(grid.best_params_) #将字符串写入文件中
        f.write("\n")   #换行 
        f.close()

