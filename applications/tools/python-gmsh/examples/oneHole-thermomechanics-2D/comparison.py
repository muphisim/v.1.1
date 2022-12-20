import pandas as pd
import matplotlib.pyplot as plt


plt.close("all")

#plt.figure()
#mydatq = pd.read_csv("muphisim/TipDisplacement_Unknown1_Mean.csv",sep=",")
#t = mydatq["Time"].values
#uy = mydatq["Value"].values
#plt.plot(t,uy,'g--',label="Muphisim")
#
#mydatq = pd.read_csv("muphisim-mpi/TipDisplacement_Unknown1_Mean_rank0.csv",sep=",")
#t = mydatq["Time"].values
#uy = mydatq["Value"].values
#plt.plot(t,uy,'k.:',label="Muphisim-mpi")
#
#plt.xlabel("Time")
#plt.ylabel("uy tip")
#plt.legend()
#plt.savefig("comparison-dispTip.pdf",)
#
#plt.figure()
#mydatq = pd.read_csv("muphisim/TipDisplacement_Unknown3_Mean.csv",sep=",")
#t = mydatq["Time"].values
#uy = mydatq["Value"].values
#plt.plot(t,uy,'g--',label="Muphisim")
#
#mydatq = pd.read_csv("muphisim-mpi/TipDisplacement_Unknown3_Mean_rank0.csv",sep=",")
#t = mydatq["Time"].values
#uy = mydatq["Value"].values
#plt.plot(t,uy,'k.:',label="Muphisim-mpi")
#
#plt.xlabel("Time")
#plt.ylabel("T tip")
#plt.legend()
#plt.savefig("comparison-T.pdf",)

plt.figure()
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
print(defo)

mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank0.csv",sep=",")
t = mydata["Time"].values
defo = mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank1.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank2.csv",sep=",")
defo = defo+ mydata["Value"].values
mydata = pd.read_csv("muphisim-mpi/All_KIN_ENERGY_Sum_rank3.csv",sep=",")
defo = defo+ mydata["Value"].values
plt.plot(t,defo,'rv',label="muphisim-mpi, kinetic")

plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig("comparison-energy.pdf",)

plt.show()
