import matplotlib.pyplot as plt
import pandas as pd

plt.close("all")
plt.figure()
d= pd.read_csv("result-fenics.csv")
x = d.values[:,0]
for cc in d.columns:
    if cc == "$t=100$" or cc == "$t=200$" or  cc == "$t=400$" or cc == "$t=10000$":
        print(cc)
        y = d[cc].values
        plt.plot(x, y, label = cc+"s, fenics")
    
p = pd.read_csv("result-100.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values-293.
plt.plot(x,y,".",label="current-100s")

p = pd.read_csv("result-200.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values-293.
plt.plot(x,y,".",label="current-200s")

p = pd.read_csv("result-400.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values-293.
plt.plot(x,y,".",label="current-400s")

p = pd.read_csv("result-10000.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values-293.
plt.plot(x,y,".",label="current-10000s")

plt.xlabel("$x$-coordinate along $y=0$")
plt.ylabel("Temperature variation $\Theta$")

plt.legend(ncol=2)

plt.savefig("comparison-with-fernics.png",dpi=1000)

plt.show()