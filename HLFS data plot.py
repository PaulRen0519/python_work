# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 12:13:49 2025

@author: admin
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.preprocessing import StandardScaler

# 加载数据
df = pd.read_csv("HLFS Full Raw Data Set.csv")

# 清洗并映射 Likert 评分
cols = ['Q5', 'Q47', 'Q22', 'Q29', 'Q14', 'Q15', 'Q43', 'Q32']
df_sub = df[cols].dropna()
likert_map = {'1': 1, '2': 2, '3': 3, '4': 4}
for col in cols:
    df_sub[col] = df_sub[col].astype(str).map(likert_map).astype(float)

# 拆分变量
X = df_sub[['Q5', 'Q22', 'Q29', 'Q14', 'Q15', 'Q43', 'Q32']]
y = df_sub['Q47']

# 标准化
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

# 建模
model = OrderedModel(y, X_scaled, distr='logit')
result = model.fit(method='bfgs')

# 可视化回归系数
coefs = result.params
conf = result.conf_int()
conf['coef'] = coefs
conf.columns = ['2.5%', '97.5%', 'coef']
conf.sort_values('coef', inplace=True)

plt.figure(figsize=(8, 5))
sns.pointplot(data=conf.reset_index(), x='coef', y='index', join=False, color='blue')
plt.axvline(0, color='gray', linestyle='--')
plt.title("Standardized Coefficients from Ordinal Logistic Regression")
plt.xlabel("Standardized Coefficient")
plt.ylabel("Variable")
plt.tight_layout()
plt.show()
