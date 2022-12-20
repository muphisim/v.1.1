import pandas as pd
import matplotlib.pyplot as plt


plt.close("all")

plt.figure()

dataFenics = pd.read_csv("fenics/dispyTip.csv",sep=",")
timeFenics = dataFenics["Time"].values
uyFenics = dataFenics["uy"].values
plt.plot(timeFenics,uyFenics,'r-',label="Fenics")

dataFenics = pd.read_csv("fenics-coarse/dispyTip.csv",sep=",")
timeFenics = dataFenics["Time"].values
uyFenics = dataFenics["uy"].values
plt.plot(timeFenics,uyFenics,'kd',label="Fenics-coarse")


mydatq = pd.read_csv("muphisim/TipDisplacement_Unknown1_Mean.csv",sep=",")
t = mydatq["Time"].values
uy = mydatq["Value"].values
plt.plot(t,uy,'g--',label="Muphisim")

mydatq = pd.read_csv("muphisim-mpi/TipDisplacement_Unknown1_Mean_rank0.csv",sep=",")
t = mydatq["Time"].values
uy = mydatq["Value"].values
plt.plot(t,uy,'k.:',label="Muphisim-mpi")

mydatq = pd.read_csv("muphisim-coarse/TipDisplacement_Unknown1_Mean.csv",sep=",")
t = mydatq["Time"].values
uy = mydatq["Value"].values
plt.plot(t,uy,'m-',label="Muphisim-coarse")

mydatq = pd.read_csv("muphisim-mpi-coarse/TipDisplacement_Unknown1_Mean_rank1.csv",sep=",")
t = mydatq["Time"].values
uy = mydatq["Value"].values
plt.plot(t,uy,'rv',label="Muphisim-mpi-coarse")


plt.xlabel("Time")
plt.ylabel("uy at the point (1,0.05,0)")
plt.legend()
plt.savefig("comparison-dispTip.pdf",)

plt.figure()

dataFenics = pd.read_csv("fenics/energy.csv",sep=",")
t = dataFenics["Time"].values
defo = dataFenics["elastic"].values
kin = dataFenics["kinetic"].values

plt.plot(t,defo,'r-',label="Fenics, defo")
plt.plot(t,kin,'g-',label="Fenics, kinetic")


mydata = pd.read_csv("muphisim/All_DEFO_ENERGY_Sum.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
plt.plot(t,defo,'k--',label="muphisim, defo")

mydata = pd.read_csv("muphisim/All_KIN_ENERGY_Sum.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
plt.plot(t,defo,'c--',label="muphisim, kinetic")

mydata = pd.read_csv("muphisim-mpi/All_DEFO_ENERGY_Sum_rank0.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_DEFO_ENERGY_Sum_rank1.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_DEFO_ENERGY_Sum_rank2.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_DEFO_ENERGY_Sum_rank3.csv",sep=",")
defo = defo+ mydata["Value"].values
plt.plot(t,defo,'m.:',label="muphisim-mpi, defo")

mydata = pd.read_csv("muphisim/All_KIN_ENERGY_Sum.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
plt.plot(t,defo,'c--',label="muphisim, kinetic")

mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank0.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank1.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank2.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank3.csv",sep=",")
defo = defo+ mydata["Value"].values
plt.plot(t,defo,'m.:',label="muphisim-mpi, kinetic")

plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig("comparison-energy.pdf",)

plt.show()
