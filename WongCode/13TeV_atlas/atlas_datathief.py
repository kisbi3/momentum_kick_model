import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit


mpl.rcParams["text.usetex"] = True

def periph(phi, a, b):
    return a*np.cos(phi)+b


def templ_276(phi, dataintegral, F, v):
    G=(dataintegral-F*norm_276)/np.pi
    return F*periph_276_data+G*(1+2*v*np.cos(2*phi))

def templ_130_G_periph(phi, G, F, v, i):
    return F*periph_130_data[i]+G

def templ_130_G_periph0(phi, G, F, v):
    return F*periph_130_data[10]+G

def templ_130_ridge_periph0(phi, G, F, v):
    return F*periph_130_data[10]+G*(1+2*v*np.cos(2*phi))

def templ_130(phi, G, F, v):
    return F*periph_130_data+G*(1+2*v*np.cos(2*phi))

Gfit=[]
Ffit = []
vfit = []
def wrap_templ_130(dataintegral):
    Gfit.clear()
    Ffit.clear()
    vfit.clear()
    def tempfunc(phi, F, v, integral = dataintegral):
        G_dist=(integral-F*norm_130)/np.pi
        Gfit.append(G_dist)
        Ffit.append(F)
        vfit.append(v)
        return templ_130(phi, G_dist, F, v)
    
    return tempfunc


fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(25, 20, forward=True)

periph_276_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[0],skiprows=1)
periph_276_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[1],skiprows=1)
periph_130_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_130_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)

plt.scatter(periph_276_phi, periph_276_data, color = 'red', s=100,label=fr'$peripheral, 2.76TeV$')
plt.scatter(periph_130_phi, periph_130_data, color = 'blue', s=100,label=fr'$peripheral, 13TeV$')

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('peripheral.png')

fig.clear()

#templ, data 등을 모두 그려보기 위해 peripheral의 0~pi 까지의 적분값을 구해보자.

#2.76TeV, 13TeV의 개수는 동일함.
periph_130_phi_avg=np.zeros(len(periph_130_phi)-1)
periph_276_phi_avg=np.zeros(len(periph_276_phi)-1)
for i in range(len(periph_130_phi)-1):
    periph_130_phi_avg[i]=(periph_130_phi[i]+periph_130_phi[i+1])/2.
    periph_276_phi_avg[i]=(periph_276_phi[i]+periph_276_phi[i+1])/2.
end1 = 2*periph_130_phi[-1]-periph_130_phi_avg[-1]
end2 = 2*periph_276_phi[-1]-periph_276_phi_avg[-1]

periph_130_phi_avg = np.append(periph_130_phi_avg, end1)
periph_276_phi_avg = np.append(periph_276_phi_avg, end2)

#적분값을 구해봅시다.
norm_276 = 0.
norm_130 = 0.
for i in range(len(periph_130_phi_avg)-1):
    if(-0.001<periph_130_phi_avg[i]<np.pi):
        norm_130 += periph_130_data[i+1]*(periph_130_phi_avg[i+1]-periph_130_phi_avg[i])
    if(-0.001<periph_276_phi_avg[i]<np.pi):
        norm_276 += periph_276_data[i+1]*(periph_276_phi_avg[i+1]-periph_276_phi_avg[i])


phi_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[0])
data_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[1])


norm = 0.

avg_13TeV_130_up=np.zeros(len(phi_13TeV_130_up)-1)
for i in range(len(phi_13TeV_130_up)-1):
    avg_13TeV_130_up[i]=(phi_13TeV_130_up[i]+phi_13TeV_130_up[i+1])/2.
end = 2*phi_13TeV_130_up[-1]-avg_13TeV_130_up[-1]

avg_13TeV_130_up = np.append(avg_13TeV_130_up, end)

#적분값을 구해봅시다.
norm_130_up = 0.
for i in range(len(avg_13TeV_130_up)-1):
    if(-0.001<avg_13TeV_130_up[i]<np.pi):
        norm_130_up += data_13TeV_130_up[i+1]*(avg_13TeV_130_up[i+1]-avg_13TeV_130_up[i])

popt_13TeV_130_up, pcov_13TeV_130_up = curve_fit(wrap_templ_130(norm_130_up), phi_13TeV_130_up, data_13TeV_130_up)

plt.scatter(phi_13TeV_130_up, data_13TeV_130_up, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_130_up))
G_periph0=np.zeros(len(phi_13TeV_130_up))
ridge_periph0=np.zeros(len(phi_13TeV_130_up))

for i in range(len(phi_13TeV_130_up)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_130_up, templ_130(phi_13TeV_130_up, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_130_up, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_130_up, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_130_up, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_130_up, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-1.55,4.7)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (130 \leq N^{rec}_{ch})$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('13TeV_130~.png')

fig.clear()



phi_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[0])
data_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[1])


norm = 0.

avg_13TeV_120_130=np.zeros(len(phi_13TeV_120_130)-1)
for i in range(len(phi_13TeV_120_130)-1):
    avg_13TeV_120_130[i]=(phi_13TeV_120_130[i]+phi_13TeV_120_130[i+1])/2.
end = 2*phi_13TeV_120_130[-1]-avg_13TeV_120_130[-1]

avg_13TeV_120_130 = np.append(avg_13TeV_120_130, end)

#적분값을 구해봅시다.
norm_120_130 = 0.
for i in range(len(avg_13TeV_120_130)-1):
    if(-0.001<avg_13TeV_120_130[i]<np.pi):
        norm_120_130 += data_13TeV_120_130[i+1]*(avg_13TeV_120_130[i+1]-avg_13TeV_120_130[i])

popt_13TeV_120_130, pcov_13TeV_120_130 = curve_fit(wrap_templ_130(norm_120_130), phi_13TeV_120_130, data_13TeV_120_130)

plt.scatter(phi_13TeV_120_130, data_13TeV_120_130, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_120_130))
G_periph0=np.zeros(len(phi_13TeV_120_130))
ridge_periph0=np.zeros(len(phi_13TeV_120_130))

for i in range(len(phi_13TeV_120_130)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_120_130, templ_130(phi_13TeV_120_130, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_120_130, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_120_130, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_120_130, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_120_130, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-1.55,4.7)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (120 \leq N^{rec}_{ch} < 130)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('13TeV_120~130.png')

fig.clear()


phi_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[0])
data_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[1])


norm = 0.

avg_13TeV_110_120=np.zeros(len(phi_13TeV_110_120)-1)
for i in range(len(phi_13TeV_110_120)-1):
    avg_13TeV_110_120[i]=(phi_13TeV_110_120[i]+phi_13TeV_110_120[i+1])/2.
end = 2*phi_13TeV_110_120[-1]-avg_13TeV_110_120[-1]

avg_13TeV_110_120 = np.append(avg_13TeV_110_120, end)

#적분값을 구해봅시다.
norm_110_120 = 0.
for i in range(len(avg_13TeV_110_120)-1):
    if(-0.001<avg_13TeV_110_120[i]<np.pi):
        norm_110_120 += data_13TeV_110_120[i+1]*(avg_13TeV_110_120[i+1]-avg_13TeV_110_120[i])

popt_13TeV_110_120, pcov_13TeV_110_120 = curve_fit(wrap_templ_130(norm_110_120), phi_13TeV_110_120, data_13TeV_110_120)

plt.scatter(phi_13TeV_110_120, data_13TeV_110_120, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_110_120))
G_periph0=np.zeros(len(phi_13TeV_110_120))
ridge_periph0=np.zeros(len(phi_13TeV_110_120))

for i in range(len(phi_13TeV_110_120)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_110_120, templ_130(phi_13TeV_110_120, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_110_120, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_110_120, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_110_120, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_110_120, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-1.55,4.7)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (110 \leq N^{rec}_{ch}<120)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('13TeV_110~120.png')

fig.clear()