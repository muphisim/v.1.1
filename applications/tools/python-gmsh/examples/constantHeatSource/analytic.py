import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.close("all")

plt.figure()

p = pd.read_csv("result.csv") 
x = p["Points_1"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="Current")

k = 10.
L = 2.
q = 1e4

T0 = 10
T1 = 50

a = (T1-T0)/L + q*L/(2*k)
b = T0

x = np.linspace(0,L,100)
y = -x**2*q/(2*k) + a*x + b

plt.plot(x,y,label="analytic")

plt.legend()
plt.show()