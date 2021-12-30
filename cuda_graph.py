import matplotlib.pyplot as plt
import numpy as np
import cupy as cp

pt=np.loadtxt('Jet.csv',delimiter=',',usecols=[0],skiprows=1)
Njetdis=np.loadtxt('Jet.csv',delimiter=',',usecols=[1],skiprows=1)

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)

pt_gpu = cp.asarray(pt)
Njetdis_gpu = cp.asarray(Njetdis)

