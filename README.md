# MethyDeep

If you just want to use our model to predict, follow the following process

```python
### load model
from keras.models import load_model
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score
new_model = load_model('finall_model.h5')

### load data
data = pd.read_csv("./yourdata")
mean = pd.read_csv("mean_30cf.csv")
std = pd.read_csv("sd_30cg.csv")
data_zcore=(data-mean)/std

#### label
data = pd.read_csv("RF_30.csv")
target_names = data.iloc[:,0].unique()
#### predict
predict = np.argmax(new_model.predict(data_zcore), axis=-1)
predict2 = []
for j in predict: 
    predict2.append(target_names[j])
predict2 ### predict result

```