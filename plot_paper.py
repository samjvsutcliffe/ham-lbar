import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate
PDF_OUTPUT = False


#plt.rc('font', family='serif', serif='Times')
## plt.rc('text', usetex=True)
#plt.rc('xtick', labelsize=8)
#plt.rc('ytick', labelsize=8)
#plt.rc('axes', labelsize=8)
plt.style.use("seaborn-paper")
plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618

def calculate_gf(disp,load):
    return integrate.trapz(load,disp)#/(0.102*0.6*13e-3)

data = pd.read_csv("load-disp.csv")


print("GF experimental:",calculate_gf(1e0*data["disp"],100e-3*data["load"]))


regex = re.compile(r'output-.*')
folders = list(filter(regex.search,os.listdir("./")))

plt.figure()
plt.plot(data["disp"],data["load"],label="FEM")
lower = pd.read_csv("lower.csv")
upper = pd.read_csv("upper.csv")
x_min = min(lower["x"].min(),upper["x"].min())
x_max = max(lower["x"].max(),upper["x"].max())
x_samples = np.linspace(x_min,x_max,100)

lower_y = np.interp(x_samples,lower["x"],1e3*lower["y"])
upper_y = np.interp(x_samples,upper["x"],1e3*upper["y"])

plt.fill_between(x_samples,lower_y,upper_y, facecolor='grey', alpha=0.5,label="Experimental")

folders = ["./output-0.5-0.0-1.0/",
           "./output-2.0-0.0-1.0/",
           "./output-4.0-0.0-1.0/"]

names = ["Coarse","Medium","Fine"]
for i,name in zip(folders,names):
    print("loading folder: ",i)
    mpm = pd.read_csv("./{}/disp.csv".format(i))
    print("GF mpm:",calculate_gf(1e3*mpm["disp"],100e-3*mpm["load"]))
    plt.plot(1e3*mpm["disp"],100e-3*mpm["load"],label=name)
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Load (N)")
plt.tight_layout()
plt.legend()
#plt.legend(["FEM","Experimental","Coarse","Medium","Fine"])
plt.savefig("paper.pdf")
plt.show()
