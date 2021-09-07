import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from io import StringIO

# mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

# pt=np.loadtxt('Jet.csv',delimiter=',',usecols=[0],skiprows=1)
jetdis=np.loadtxt('Jet_pTdis.csv',delimiter=',',usecols=[1],skiprows=1)
pt=np.loadtxt('Jet_pTdis.csv',delimiter=',',usecols=[0],skiprows=1)
Ridgedis=np.loadtxt('Ridge_pTdis.csv',delimiter=',',usecols=[1],skiprows=1)

CMS7TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[0],skiprows=13,max_rows=34)
CMS7TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[3],skiprows=13,max_rows=34)
CMS7TeVptdis_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[4],skiprows=13,max_rows=34)
CMS7TeVptdis_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[5],skiprows=13,max_rows=34)

# print(CMS7TeVpt)

plt.errorbar(CMS7TeVpt,CMS7TeVptdis, yerr=(CMS7TeVptdis_error1,abs(CMS7TeVptdis_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
plt.plot(pt,jetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
plt.plot(pt,Ridgedis+0.632*jetdis, color="black",linewidth=7,label=r'$p+p,Jet+Ridge$')
plt.plot(pt,Ridgedis, color="red",linewidth=4,linestyle = '-',label=r'$p+p,Ridge$')


plt.minorticks_on()
plt.yscale('log')
ax.axis([0,4,0.0004,20])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('pTdis.png')
# fig.savefig('Initial_Parton.png')
# plt.show()