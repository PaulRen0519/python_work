import pandas as pd

df = pd.read_excel('HLFS Full Raw Data Set.xlsx')
df = df.dropna()
x = df['How often do you have breakfast during the academic semester?']
y = df['On average, how many hours per day do you use your cellphone?']

import statsmodels.api as sm
x = sm.add_constant(x)
model = sm.OLS(y, x).fit()
print(model.summary())