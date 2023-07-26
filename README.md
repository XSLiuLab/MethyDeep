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
  * [validata_dataset](./code/validata_dataset) : Performance of MethyDeep in validata dataset
  * [model_comparsion](./code/model_comparsion/) : Performance of MethyDeep and other  models on validata dataset
  * [plot](./code/plot) : All drawing codes in this study
* [data](./data) : [RF_feature45.csv](./data/RF_feature45.csv) : The features of MethyDeep ; [targer_names.csv](https://github.com/XSLiuLab/MethyDeep/blob/main/data/targer_names.csv) : model lables
* [MethyDeep](./code/MethyDeep/) :  Finalized model

## Brief use

If you just want to use our model to predict, follow the following process

~~~python
### load model
from keras.models import load_model
import pandas as pd
import numpy as np
new_model = load_model('./MethyDeep')
feature = pd.read_csv("./data/RF_feature45.csv")
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
Ning W, Wu T, Wu C, Wang S, Tao Z, Wang G, Zhao X, Diao K, Wang J, Chen J, Chen F, Liu XS. Accurate prediction of pan-cancer types using machine learning with minimal number of DNA methylation sites. J Mol Cell Biol. 2023 Apr 10:mjad023. doi: 10.1093/jmcb/mjad023. Epub ahead of print. PMID: 37037781.

***

**Cancer Biology Group @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**
