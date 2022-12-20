import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams.update({'font.size': 12})


def getData(srcDir, destDir):
    data = []
    columns = None
    tests = ["plate-0","plate-1","plate-2","plate-3"]
    for te in tests:
        fname = os.path.join(srcDir,te,"errorData_step1.csv")
        if os.path.exists(fname):
            p = pd.read_csv(fname,sep=",")
            if columns is None:
                columns = p.columns
            data.append(p.values.tolist()[0])
            print(data)
    p = pd.DataFrame(data,columns=columns)
    p.to_csv(os.path.join(destDir,"allResults.csv"),index=False)
    return p




plt.close("all")

f, (ax1, ax2) = plt.subplots(1, 2, figsize=[6.4*2, 4.8] )


p = getData("FEM","FEM")
ax1.plot(p.values[:,0]**0.5,p.values[:,2],"r.-",markersize=8,label="FEM")
ax2.plot(p.values[:,0]**0.5,p.values[:,3],"r.-",markersize=8,label="FEM")


p = getData("MM","MM")
ax1.plot(p.values[:,0]**0.5,p.values[:,2],"gs--",markersize=8,label="MM")
ax2.plot(p.values[:,0]**0.5,p.values[:,3],"gs--",markersize=8,label="MM")

ax2.legend(ncol=1, loc=(1.01, 0))
ax1.set_xscale("log")
ax1.set_yscale("log")

ax2.set_xscale("log")
ax2.set_yscale("log")

ax1.set_title(r"(A) $L_2$-norm error",fontsize=16)
ax1.set_xlabel(r"$\sqrt{N_{el}}$",fontsize=16)
ax1.set_ylabel(r"$\varepsilon_{L_2}$",fontsize=16)

ax2.set_title(r"(B) $H_1$-norm error",fontsize=16)
ax2.set_xlabel(r"$\sqrt{N_{el}}$",fontsize=16)
ax2.set_ylabel(r"$\varepsilon_{H_1}$",fontsize=16)

plt.savefig("ErrorPlate.eps",bbox_inches='tight',dpi=1000)
os.system("epstopdf ErrorPlate.eps")



plt.show()


