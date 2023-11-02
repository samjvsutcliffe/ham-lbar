import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.figure()
print("plotting")
df = pd.read_csv("output/terminus_position.csv")
time = df["Time (s)"].values
position = df["Terminus position"].values
time_days = time/(24*60*60)

stokes_data = pd.read_csv("stokes.csv")


plt.plot(time_days[1:],position[1:] - position[1],label="MPM")
plt.plot(stokes_data["time"] * 30, stokes_data["disp"],label="Stokes")
plt.xlabel("Time (d)")
plt.legend()

plt.figure()
plt.plot(time_days,np.gradient(position,time_days),label="MPM")
time_days = stokes_data["time"] * 30
position = stokes_data["disp"]
plt.plot(time_days,np.gradient(position,time_days),label="Stokes")
plt.legend()
#for f in [6,7]:
#    print("plotting {}".format(f));
#    df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#    plt.plot(df["Time (s)"],df["Terminus position"]*(10**f),label=f)
#f = 6
#df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#df["Terminus position"] -= df["Terminus position"][0]
#plt.plot(df["Time (s)"],df["Terminus position"],label=f)
#f = 7
#df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#df["Terminus position"] -= df["Terminus position"][0]
#plt.plot(df["Time (s)"],df["Terminus position"]*10,label=f)
plt.show()

