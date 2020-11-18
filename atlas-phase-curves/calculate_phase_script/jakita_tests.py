import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

df=pd.read_csv("jakita_tests.csv")
print(df)


fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

ax1.scatter(df['Threads'],df['Objects/s'])

ax1.set_xlabel("Threads")
ax1.set_ylabel("Objects/s")

plt.show()
