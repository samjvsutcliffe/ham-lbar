import matplotlib.pyplot as plt
import pandas as pd

from scipy import integrate
def calculate_gf(disp,load):
    return integrate.trapz(load,disp)#/(0.102*0.6*13e-3)

data = pd.read_csv("load-disp.csv")


print("GF experimental:",calculate_gf(1e0*data["disp"],100e-3*data["load"]))

mpm = pd.read_csv("output/disp.csv")

plt.plot(data["disp"],data["load"],label="Data")
print("GF mpm:",calculate_gf(1e3*mpm["disp"],100e-3*mpm["load"]))
plt.plot(1e3*mpm["disp"], 100e-3*mpm["load"],label="MPM")
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend()
plt.show()
