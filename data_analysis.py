import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.preprocessing import StandardScaler

# 1. 读取数据
df = pd.read_csv("HLFS Full Raw Data Set.csv")

# 2. 选择变量并去除缺失值
cols = ['Q5', 'Q47', 'Q22', 'Q29', 'Q14', 'Q15', 'Q43', 'Q32']
df_sub = df[cols].dropna()

# 3. 编码 Likert 量表（按你问卷中的顺序）
likert_map = {
    '1': 1,
    '2': 2,
    '3': 3,
    '4': 4,
    }

for col in cols:
    df_sub[col] = df_sub[col].map(likert_map).fillna(df_sub[col])

# 4. 建模数据
X = df_sub[['Q5', 'Q22', 'Q29', 'Q14', 'Q15', 'Q43', 'Q32']]
y = df_sub['Q47']

# 5. 标准化自变量
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

# 6. 构建有序Logistic回归模型
model = OrderedModel(y, X_scaled, distr='logit')
result = model.fit(method='bfgs')

# 7. 打印回归结果
print(result.summary())

# 8. 可视化：系数图
coefs = result.params
conf = result.conf_int()
conf['coef'] = coefs
conf.columns = ['2.5%', '97.5%', 'coef']
conf.sort_values('coef', inplace=True)

plt.figure(figsize=(8,5))
sns.pointplot(data=conf.reset_index(), x='coef', y='index', join=False)
plt.axvline(0, color='gray', linestyle='--')
plt.title("变量对抑郁频率的影响（标准化回归系数）")
plt.xlabel("标准化回归系数")
plt.ylabel("变量")
plt.tight_layout()
plt.show()

