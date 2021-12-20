import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit


mpl.rcParams["text.usetex"] = True

errors = 0

f = open('deviation.csv','w',newline='')
wr=csv.writer(f)

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)
# fig.set_size_inches(50, 30, forward=True)

PbPb276_pttrig03_1_deltaphi=np.loadtxt('long_pttrig03-1.csv',delimiter=',',usecols=[0])
PbPb276_pttrig03_1=np.loadtxt('long_pttrig03-1.csv',delimiter=',',usecols=[1])
pttrig03_1min=min(PbPb276_pttrig03_1)
PbPb276_pttrig1_2_deltaphi=np.loadtxt('long_pttrig1-2.csv',delimiter=',',usecols=[0])
PbPb276_pttrig1_2=np.loadtxt('long_pttrig1-2.csv',delimiter=',',usecols=[1])
pttrig1_2min=min(PbPb276_pttrig1_2)
PbPb276_pttrig2_4_deltaphi=np.loadtxt('long_pttrig2-4.csv',delimiter=',',usecols=[0])
PbPb276_pttrig2_4=np.loadtxt('long_pttrig2-4.csv',delimiter=',',usecols=[1])
pttrig2_4min=min(PbPb276_pttrig2_4)
PbPb276_pttrig4_6_deltaphi=np.loadtxt('long_pttrig4-6.csv',delimiter=',',usecols=[0])
PbPb276_pttrig4_6=np.loadtxt('long_pttrig4-6.csv',delimiter=',',usecols=[1])
pttrig4_6min=min(PbPb276_pttrig4_6)
PbPb276_pttrig6_12_deltaphi=np.loadtxt('long_pttrig6-12.csv',delimiter=',',usecols=[0])
PbPb276_pttrig6_12=np.loadtxt('long_pttrig6-12.csv',delimiter=',',usecols=[1])
pttrig6_12min=min(PbPb276_pttrig6_12)

resultphi_031=np.loadtxt('phiCorrelation_pt0-1.csv',delimiter=',',usecols=[0],skiprows=1)
pttrig_031=np.loadtxt('phiCorrelation_pt0-1.csv',delimiter=',',usecols=[2],skiprows=1)
phicorrmin031 = min(pttrig_031)
resultphi_12=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1)
pttrig_12=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1)
phicorrmin12 = min(pttrig_12)
resultphi_24=np.loadtxt('phiCorrelation_pt2-4.csv',delimiter=',',usecols=[0],skiprows=1)
pttrig_24=np.loadtxt('phiCorrelation_pt2-4.csv',delimiter=',',usecols=[2],skiprows=1)
phicorrmin24 = min(pttrig_24)
resultphi_46=np.loadtxt('phiCorrelation_pt4-6.csv',delimiter=',',usecols=[0],skiprows=1)
pttrig_46=np.loadtxt('phiCorrelation_pt4-6.csv',delimiter=',',usecols=[2],skiprows=1)
phicorrmin46 = min(pttrig_46)
resultphi_612=np.loadtxt('phiCorrelation_pt6-12.csv',delimiter=',',usecols=[0],skiprows=1)
pttrig_612=np.loadtxt('phiCorrelation_pt6-12.csv',delimiter=',',usecols=[2],skiprows=1)
phicorrmin612 = min(pttrig_612)

pttrig031_czyam = PbPb276_pttrig03_1-pttrig03_1min
pttrig12_czyam = PbPb276_pttrig1_2-pttrig1_2min
pttrig24_czyam = PbPb276_pttrig2_4-pttrig2_4min
pttrig46_czyam = PbPb276_pttrig4_6-pttrig4_6min
pttrig612_czyam = PbPb276_pttrig6_12-pttrig6_12min

pttrig_031 = pttrig_031 - phicorrmin031
pttrig_12 = pttrig_12 - phicorrmin12
pttrig_24 = pttrig_24 - phicorrmin24
pttrig_46 = pttrig_46 -phicorrmin46
pttrig_612 = pttrig_612 - phicorrmin612

plt.scatter(PbPb276_pttrig03_1_deltaphi, pttrig031_czyam, color = "black", s = 80)


plt.plot(resultphi_031, pttrig_031, color = 'blue', linewidth=7, linestyle = '-')

plt.xlabel(r'$\Delta\phi$',size=70)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

plt.title(r'$0.3<p^{trig}_T<1$', fontsize = 75)
plt.minorticks_on()

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_long_pttrig03-1.png')


fig.clear()


plt.scatter(PbPb276_pttrig1_2_deltaphi, pttrig12_czyam, color = "black", s = 80)
plt.plot(resultphi_12, pttrig_12, color = 'blue', linewidth=7, linestyle = '-')

plt.xlabel(r'$\Delta\phi$',size=70)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

plt.title(r'$1<p_T^{trig}<2$', fontsize = 75)
plt.minorticks_on()

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()


fig.savefig('1D_PhiCorrelation_long_pttrig1-2.png')

fig.clear()

plt.scatter(PbPb276_pttrig2_4_deltaphi, pttrig24_czyam, color = "black", s = 80)
plt.plot(resultphi_24, pttrig_24, color = 'blue', linewidth=7, linestyle = '-')

plt.xlabel(r'$\Delta\phi$',size=70)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

plt.title(r'$2<p_T^{trig}<4$', fontsize = 75)
plt.minorticks_on()

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()


fig.savefig('1D_PhiCorrelation_long_pttrig2-4.png')

fig.clear()


plt.scatter(PbPb276_pttrig4_6_deltaphi, pttrig46_czyam, color = "black", s = 80)
plt.plot(resultphi_46, pttrig_46, color = 'blue', linewidth=7, linestyle = '-')


plt.xlabel(r'$\Delta\phi$',size=70)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

plt.title(r'$4<p_T^{trig}<6$', fontsize = 75)
plt.minorticks_on()

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_long_pttrig4-6.png')

fig.clear()

plt.scatter(PbPb276_pttrig6_12_deltaphi, pttrig612_czyam, color = "black", s = 80)
plt.plot(resultphi_612, pttrig_612, color = 'blue', linewidth=7, linestyle = '-')

plt.xlabel(r'$\Delta\phi$',size=70)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

plt.title(r'$6<p_T^{trig}<12$', fontsize = 75)
plt.minorticks_on()

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_long_pttrig6-12.png')

fig.clear()


