import pandas as pd
import matplotlib.pyplot as plt


plt.figure()

a = pd.read_csv('perpFlap-explicit/solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
plt.plot(a['Time'].values, a['Displacement0'].values,'g-',label='MuPhiSim-explicit')

a = pd.read_csv('perpFlap-explicit-mpi/solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
plt.plot(a['Time'].values, a['Displacement0'].values,'m--',label='MuPhiSim-explicit-parallel')


a = pd.read_csv('perpFlap-mpi/solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
plt.plot(a['Time'].values, a['Displacement0'].values,'r-',label='MuPhiSim-parallel')

a = pd.read_csv('perpFlap/solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
plt.plot(a['Time'].values, a['Displacement0'].values,'g:',label='MuPhiSim-sequential')

#a = pd.read_csv('perpFlap-smallStep/solid-MuPhiSim/precice-Solid-watchpoint-Flap-Tip.log',delim_whitespace=True)
#plt.plot(a['Time'].values, a['Displacement0'].values,'c:',label='MuPhiSim-sequential - small step')


plt.legend()

plt.xlabel("Time")
plt.ylabel("Displacement0")

plt.savefig("comparison.pdf")

plt.show()
