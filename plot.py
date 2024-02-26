import matplotlib.pyplot as plt
import pandas as pd
import os
import re

from scipy import integrate
def calculate_gf(disp,load):
    return integrate.trapz(load,disp)#/(0.102*0.6*13e-3)

data = pd.read_csv("load-disp.csv")


print("GF experimental:",calculate_gf(1e0*data["disp"],100e-3*data["load"]))


regex = re.compile(r'output-.*')
folders = list(filter(regex.search,os.listdir("./")))

for i in folders:
    print("loading folder: ",i)
    plt.figure()
    mpm = pd.read_csv("./{}/disp.csv".format(i))
    plt.plot(data["disp"],data["load"],label="Data")
    print("GF mpm:",calculate_gf(1e3*mpm["disp"],100e-3*mpm["load"]))
    plt.plot(1e3*mpm["disp"],100e-3*mpm["load"],label="MPM")
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Load (N)")
    plt.legend()
plt.show()
