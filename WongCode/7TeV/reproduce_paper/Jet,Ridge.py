import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)
# fig.set_size_inches(20, 25, forward=True)

pt=np.loadtxt('7TeV_Jet.csv',delimiter=',',usecols=[0],skiprows=1)
Njetdis=np.loadtxt('7TeV_Jet.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge=np.loadtxt('7TeV_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_remove_Ei=np.loadtxt('7TeV_Ridge_remove_Ei.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_remove_allE=np.loadtxt('7TeV_Ridge_remove_allE.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_asdf=np.loadtxt('Ridge_Jet_asdf.csv',delimiter=',',usecols=[1],skiprows=1)
# piRidge=np.loadtxt('md=mpi_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)

Jetpt_Paper=np.loadtxt('7TeV_Jet_Paper.csv',delimiter=',',usecols=[0])
Jet_Paper=np.loadtxt('7TeV_Jet_Paper.csv',delimiter=',',usecols=[1])
Ridgept_Paper=np.loadtxt('7TeV_Ridge_Paper.csv',delimiter=',',usecols=[0])
Ridge_Paper=np.loadtxt('7TeV_Ridge_Paper.csv',delimiter=',',usecols=[1])

hepdata_7tevpt=np.loadtxt('7TeV_HepData.csv',delimiter=',',usecols=[0],skiprows=1)
hepdata_7tev=np.loadtxt('7TeV_HepData.csv',delimiter=',',usecols=[1],skiprows=1)
hepdata_7tev_error1=tuple(np.loadtxt('7TeV_HepData.csv',delimiter=',',usecols=[2],skiprows=1))
hepdata_7tev_error2=tuple(abs(np.loadtxt('7TeV_HepData.csv',delimiter=',',usecols=[3],skiprows=1)))
err=[hepdata_7tev_error1,hepdata_7tev_error2]
# print(err)

# for i in range(len):
#     a=hepdata_7tev_error1[i]
#     b=hepdata_7tev_error2[i]
#     print(i)
#     err.append((a),(b))
# print(err)
# RefRidgempipt=np.loadtxt('Au+Au,Ridge,md=mpi.csv',delimiter=',',usecols=[0])
# RefRidgempi=np.loadtxt('Au+Au,Ridge,md=mpi.csv',delimiter=',',usecols=[1])
# RefRidgept=np.loadtxt('Au+Au,Ridge.csv',delimiter=',',usecols=[0])
# RefRidge=np.loadtxt('Au+Au,Ridge.csv',delimiter=',',usecols=[1])
# ptJeongseok=np.loadtxt('JeongSeok/Ridge_Jet.csv',delimiter=',',usecols=[0],skiprows=1)
# RidgeJetJeongseok=np.loadtxt('JeongSeok/Ridge_Jet.csv',delimiter=',',usecols=[1],skiprows=1)

# Ridge_remove_Ei = np.loadtxt('Ridge_remove_Ei.csv',delimiter=',',usecols=[1],skiprows=1)

plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
plt.plot(pt,Njetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
# plt.plot(pt,(Ridge*8/3)+0.632*Njetdis, color="black",linewidth=7,label=r'$pp,Jet+Ridge$')
plt.plot(pt,Ridge, color="red",linewidth=7,linestyle = '--',label=r'$p+p,Ridge$')
# plt.plot(pt,Ridge_remove_Ei, color="green",linewidth=7,linestyle = '--',label=r'$p+p,Ridge remove E_i$')
# plt.plot(pt,Ridge_remove_allE, color="blue",linewidth=7,linestyle = '--',label=r'$p+p,Ridge remove E/E_i$')
# plt.plot(pt,Ridge_asdf, color="black",linewidth=7,linestyle = '--',label=r'$p+p,Ridge remove E/E_i$')
# plt.plot(pt,Ridge_remove_Ei, color="red",linewidth=7,linestyle = '--',label=r'$Au+Au,Ridge remove$')

# plt.plot(ptJeongseok,RidgeJetJeongseok, color="blue",linewidth=7,linestyle = '--',label=r'$Au+Au,Ridge,Jeongseok$')

# plt.plot(pt,piRidge*8/3, color="red",linewidth=7,linestyle = ':', label=r'$Au+Au,Ridge, m_{d}=m_{\pi}$')
# plt.errorbar(hepdata_7tevpt,hepdata_7tev, yerr=err,fmt='s',linewidth=7,label=r'pp,7TeV CMS')
plt.scatter(Jetpt_Paper,Jet_Paper, color="blue",s=200, label=r'$p+p,Jet[Paper]$')
plt.scatter(Ridgept_Paper,Ridge_Paper, color="red",s=200, label=r'$p+p,Ridge[Paper]$')

plt.minorticks_on()
plt.yscale('log')
ax.axis([0,4,0.0004,12])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('7TeV_Jet,Ridge.png')
# plt.show()