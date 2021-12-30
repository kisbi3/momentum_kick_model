import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

# date = np.loadtxt('random_number_data.csv', delimiter=',', usecols=[0], skiprows=1)

hour = np.loadtxt('random_number_data.csv', delimiter=',', usecols=[1], skiprows=1)
minute = np.loadtxt('random_number_data.csv', delimiter=',', usecols=[2], skiprows=1)
second = np.loadtxt('random_number_data.csv', delimiter=',', usecols=[3], skiprows=1)
date = np.arange(len(hour))

hour = hour-12
second = second/100
minute = minute+second
minute = minute/100
hour = hour+minute

plt.plot(date, hour, linewidth = 7)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')
plt.minorticks_on()
plt.grid(color='silver',linestyle=':',linewidth=5)

plt.tight_layout

fig.savefig('random.png')

fig.clear()