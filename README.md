# **Accurate prediction of pan-cancer types using machine learning with minimal number of DNA methylation sites**

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![version](https://img.shields.io/badge/version-dev-green.svg)](https://shields.io/)

<details>
<summary>Table of content</summary>

## Table of content
   * [Overview](#Overview)
   * [Contents](#Contents)
   * [Citation](#Citation)
   * [Acknowledgement](#Acknowledgement)
   * [LICENSE](#License)

</details>

----

## Overview

This repository provides the analysis reports, code and data for readers who are interest in this project and make it easier to reproduce the whole analysis procedure.

## Contents

* [code](./code) 
  * [data_process](./code/data_process/) : Preprocessing of original data and generation of training data
  * [Difference analysis](./code/Difference-analysis/) : Difference analysis between champ and wilcox test
  * [ML](./code/ML/) : Comparison of six machine learning methods
  * [DNN_Gridsearch](./code/DNN_Gridsearch/) : Neural network for super parameter selection
  * [DNN](./code/DNN/) : Training of Neural Network
  * [metastasis_dataset](./code/metastasis_dataset) : Performance of neural network and machine learning in metastatic cancer dataset
  * [model_comparsion](./code/model_comparsion/) : Performance of methydeep and other cancer models on primary cancer dataset
* [data](./data) : [Xgboost_feature_11_23.csv](https://github.com/XSLiuLab/MethyDeep/blob/main/data/Xgboost_feature_11_23.csv) : Data features ; [targer_names.csv](https://github.com/XSLiuLab/MethyDeep/blob/main/data/targer_names.csv) : model lables
* [methydeep](./code/methydeep/) :  Finalized model

## Brief use

If you just want to use our model to predict, follow the following process

~~~python
### load model
from keras.models import load_model
import pandas as pd
import numpy as np
new_model = load_model('./code/methydeep')
feature = pd.read_csv("./data/XGboost_feature_11_23.csv")
target_names = pd.read_csv("./data/targer_names.csv")
### load data
data = pd.read_csv("./yourdata")
data = data[feature.iloc[0:30,0]]
#### predict
predict = np.argmax(new_model.predict(data), axis=-1)
predict2 = []
for j in predict: 
    predict2.append(target_names[j])
predict2 ### predict result
~~~

## Citation



## Acknowledgement

We thank ShanghaiTech University High Performance Computing Public Service Platform for computing services.This work was supported by Shanghai Science and Technology Commission (21ZR1442400), the National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.

## License

***

**Cancer Biology Group @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**