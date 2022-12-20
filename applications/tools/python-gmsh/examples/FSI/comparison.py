import pandas as pd
import matplotlib.pyplot as plt

plt.close("all")

plt.figure(figsize=(15,10))

#a = pd.read_csv('solid-MuPhiSim-qs/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
#plt.plot(a['Time'].values, a['Displacement0'].values,'g.--',label='MuPhiSim-implicit')

a = pd.read_csv('solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
plt.plot(a['Time'].values, a['Displacement0'].values,'b+',label='MuPhiSim, small densitiy')


plt.legend(ncol=2, )
plt.xlabel("Time")
plt.ylabel("Tip Displacement x")
plt.savefig("comparison.pdf", bbox_inches='tight')

plt.show()

