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
Njetdis=np.loadtxt('Jet.csv',delimiter=',',usecols=[1],skiprows=1)
pt=np.loadtxt('Ridge2_phi.csv',delimiter=',',usecols=[0],skiprows=1)
Ridge=np.loadtxt('Ridge2_phi.csv',delimiter=',',usecols=[1],skiprows=1)
init_Ridge=np.loadtxt('Initial_parton.csv',delimiter=',',usecols=[1],skiprows=1)

Ohno1=np.loadtxt('Ridge_2.csv',delimiter=',',usecols=[0],skiprows=1)
Ohno2=np.loadtxt('Ridge_2.csv',delimiter=',',usecols=[1],skiprows=1)

# piRidge=np.loadtxt('md=mpi_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)
RefJetRidgept=np.loadtxt('AuAuJet+Ridge.csv',delimiter=',',usecols=[0])
RefJetRidge=np.loadtxt('AuAuJet+Ridge.csv',delimiter=',',usecols=[1])
# RefRidgempipt=np.loadtxt('Au+Au,Ridge,md=mpi.csv',delimiter=',',usecols=[0])
# RefRidgempi=np.loadtxt('Au+Au,Ridge,md=mpi.csv',delimiter=',',usecols=[1])
RefRidgept=np.loadtxt('Au+Au,Ridge.csv',delimiter=',',usecols=[0])
RefRidge=np.loadtxt('Au+Au,Ridge.csv',delimiter=',',usecols=[1])
# ptJeongseok=np.loadtxt('JeongSeok/Ridge_Jet.csv',delimiter=',',usecols=[0],skiprows=1)
# RidgeJetJeongseok=np.loadtxt('JeongSeok/Ridge_Jet.csv',delimiter=',',usecols=[1],skiprows=1)



# Ridge_remove_Ei = np.loadtxt('Ridge_remove_Ei.csv',delimiter=',',usecols=[1],skiprows=1)

plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
plt.plot(pt,Njetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
plt.plot(pt,Ridge+0.632*Njetdis, color="black",linewidth=7,label=r'$Au+Au,Jet+Ridge$')
plt.plot(pt,Ridge, color="red",linewidth=4,linestyle = '-',label=r'$Au+Au,Ridge$')
# plt.plot(pt,init_Ridge, color="black",linewidth=7,linestyle = '--',label=r'$Au+Au,Ridge(initial)$')
# plt.plot(pt,Ridge_remove_Ei, color="red",linewidth=7,linestyle = '--',label=r'$Au+Au,Ridge remove$')

# plt.plot(ptJeongseok,RidgeJetJeongseok, color="blue",linewidth=7,linestyle = '--',label=r'$Au+Au,Ridge,Jeongseok$')

# plt.plot(pt,piRidge*8/3, color="red",linewidth=7,linestyle = ':', label=r'$Au+Au,Ridge, m_{d}=m_{\pi}$')
# plt.scatter(RefJetRidgept,RefJetRidge, color="green",s=200, label=r'$Au+Au,Jet+Ridge[Paper]$')
# plt.scatter(RefRidgempipt,RefRidgempi, color="black",s=200, label=r'$Au+Au,Ridge,md=mpi[Paper]$')
# plt.scatter(RefRidgept,RefRidge, color="deeppink",s=200, label=r'$Au+Au,Ridge[Paper]$')

Wongpt_ridge=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[2],skiprows=2, max_rows=41)
WongdNptdpt_ridge=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[4],skiprows=2, max_rows=41)

print(Wongpt_ridge)
print(WongdNptdpt_ridge)

plt.scatter(Wongpt_ridge,WongdNptdpt_ridge, color="red",s=200, label=r'$Au+Au,Ridge[Wong]$')


Wongpt_jet=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[2],skiprows=88, max_rows=41)
WongdNptdpt_jet=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[4],skiprows=88, max_rows=41)

print(Wongpt_jet)
print(WongdNptdpt_jet)

plt.scatter(Wongpt_jet,WongdNptdpt_jet, color="blue",s=200, label=r'$p+p,jet[Wong]$')


Wongpt_ridge_jet=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[2],skiprows=131, max_rows=41)
WongdNptdpt_ridge_jet=np.loadtxt('pkdnptdptden.txt',delimiter=' ',usecols=[4],skiprows=131, max_rows=41)

print(Wongpt_ridge_jet)
print(WongdNptdpt_ridge_jet)

# plt.scatter(Wongpt_ridge_jet,WongdNptdpt_ridge_jet, color="black",s=200, label=r'$Au+Au,Ridge+jet[Wong]$')
plt.scatter(Wongpt_ridge_jet,WongdNptdpt_ridge+0.632*WongdNptdpt_jet, color="black",s=200, label=r'$Au+Au,Ridge+jet[Wong]$')

# plt.scatter(Ohno1,Ohno2, color="green",s=100, label=r'$asdf$')

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

fig.savefig('Fig2_phi.png')
# fig.savefig('Initial_Parton.png')
# plt.show()