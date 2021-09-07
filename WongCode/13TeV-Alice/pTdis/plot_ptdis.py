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
alicejetdis=np.loadtxt('pTdis_alice.csv',delimiter=',',usecols=[2],skiprows=1)
alicept=np.loadtxt('pTdis_alice.csv',delimiter=',',usecols=[0],skiprows=1)
aliceRidgedis=np.loadtxt('pTdis_alice.csv',delimiter=',',usecols=[1],skiprows=1)

cmsjetdis=np.loadtxt('pTdis_cms.csv',delimiter=',',usecols=[2],skiprows=1)
cmspt=np.loadtxt('pTdis_cms.csv',delimiter=',',usecols=[0],skiprows=1)
cmsRidgedis=np.loadtxt('pTdis_cms.csv',delimiter=',',usecols=[1],skiprows=1)

Alice13TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=5)
Alice13TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=5)
Alice13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=5)
Alice13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=5)
Alice13TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[6],skiprows=12,max_rows=5)
Alice13TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[7],skiprows=12,max_rows=5)

Alice13TeVptdis_error1 = pow(Alice13TeVptdis_error1_stat*Alice13TeVptdis_error1_stat+Alice13TeVptdis_error1_sys*Alice13TeVptdis_error1_sys,0.5)
Alice13TeVptdis_error2 = pow(Alice13TeVptdis_error2_stat*Alice13TeVptdis_error2_stat+Alice13TeVptdis_error2_sys*Alice13TeVptdis_error2_sys,0.5)

# print(Alice13TeVptdis_error1)
# print(Alice13TeVptdis_error2)

CMS13TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[0],skiprows=14,max_rows=9)
CMS13TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[1],skiprows=14,max_rows=9)
CMS13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=9)
CMS13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=9)
CMS13TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=9)
CMS13TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=9)

CMS13TeVptdis_error1 = pow(CMS13TeVptdis_error1_stat*CMS13TeVptdis_error1_stat+CMS13TeVptdis_error1_sys*CMS13TeVptdis_error1_sys,0.5)
CMS13TeVptdis_error2 = pow(CMS13TeVptdis_error2_stat*CMS13TeVptdis_error2_stat+CMS13TeVptdis_error2_sys*CMS13TeVptdis_error2_sys,0.5)

# yerr=(Alice13TeVptdis_error1,abs(Alice13TeVptdis_error2))
plt.errorbar(Alice13TeVpt, Alice13TeVptdis, yerr=(abs(Alice13TeVptdis_error2),Alice13TeVptdis_error1), color="blue",markersize=17,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV, \, ALICE, \,1.6<\vert\Delta\eta\vert<1.8$',capsize=10)
plt.errorbar(CMS13TeVpt, CMS13TeVptdis, yerr=(abs(CMS13TeVptdis_error2),CMS13TeVptdis_error1), color="black",markersize=17,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV, \, CMS, \,2.0<\vert\Delta\eta\vert<4.0$',capsize=10)
# plt.plot(alicept,aliceRidgedis, color="blue",linewidth=5,linestyle = '-',label=r'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8$')
plt.scatter(alicept,aliceRidgedis, color="blue",s=1000,linewidths=4,facecolors='none',label=r'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8$')
# plt.plot(alicept,alicejetdis, color="blue",linewidth=5,linestyle = '--',label=r'$pp,Jet,1.6<\vert\Delta\eta\vert<1.8$')
# plt.plot(cmspt,cmsRidgedis, color="black",linewidth=5,linestyle = '-',label=r'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0$')
plt.scatter(cmspt,cmsRidgedis, color="black",s=1000,linewidths=4,facecolors='none',label=r'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0$')

plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
# plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
plt.ylabel(r'$Y^{Ridge}$',size=50)
# plt.plot(pt,jetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
# plt.plot(pt,Ridgedis+0.632*jetdis, color="black",linewidth=7,label=r'$p+p,Jet+Ridge$')

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([0,4,0.0004,20])

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

# fig.clear()

# Alice13TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=5)
# Alice13TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=5)
# Alice13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=5)
# Alice13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=5)
# Alice13TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[6],skiprows=12,max_rows=5)
# Alice13TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[7],skiprows=12,max_rows=5)