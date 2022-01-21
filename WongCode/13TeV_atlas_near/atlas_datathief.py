import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

ridge_parameters = open('./result/ridge_parameters.csv','w',newline='')
wr = csv.writer(ridge_parameters)  #계수 저장
wr.writerow(['Multiplicity','F','G','v_22'])


mpl.rcParams["text.usetex"] = True

def periph(phi, a, b, c, d):
    return a*np.cos(phi)+b*np.cos(2*phi)+c*np.cos(3*phi)+d


def templ_276(phi, dataintegral, F, v):
    G=(dataintegral-F*norm_276)/np.pi
    return F*periph_276_data+G*(1+2*v*np.cos(2*phi))

def templ_130_G_periph(phi, G, F, v, i):
    return F*periph_130_data[i]+G

def templ_130_G_periph0(phi, G, F, v):
    return F*periph_130_data[6]+G

def templ_130_ridge_periph0(phi, G, F, v):
    return F*periph_130_data[6]+G*(1+2*v*np.cos(2*phi))

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
        # G_dist=(integral-F*norm_130)*2/np.pi
        G_dist=(integral-F*norm_130)*(1/0.96)
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

popt_periph_130TeV, pcov_periph_130TeV = curve_fit(periph, periph_130_phi, periph_130_data)
popt_periph_276TeV, pcov_periph_276TeV = curve_fit(periph, periph_276_phi, periph_276_data)

plt.scatter(periph_276_phi, periph_276_data, color = 'red', s=100,label=fr'$peripheral, 2.76TeV$')
plt.scatter(periph_130_phi, periph_130_data, color = 'blue', s=100,label=fr'$peripheral, 13TeV$')

plt.plot(periph_130_phi, periph(periph_130_phi, *popt_periph_130TeV), color = 'blue', linewidth=5, linestyle = '-')
plt.plot(periph_276_phi, periph(periph_276_phi, *popt_periph_276TeV), color = 'red', linewidth=5, linestyle = '-')

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/peripheral.png')

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

#적분값을 구해봅시다.(peripheral integrate)
norm_276 = 0.
norm_130 = 0.
for i in range(len(periph_130_phi_avg)-1):
    if(-0.001<periph_130_phi_avg[i]):
        norm_130 += periph_130_data[i+1]*(periph_130_phi_avg[i+1]-periph_130_phi_avg[i])
        # print(periph_130_phi_avg[i+1], periph_130_phi_avg[i], periph_130_data[i+1])
    if(-0.001<periph_276_phi_avg[i]):
        norm_276 += periph_276_data[i+1]*(periph_276_phi_avg[i+1]-periph_276_phi_avg[i])


phi_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[0])
data_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[1])


norm = 0.

avg_13TeV_90_up=np.zeros(len(phi_13TeV_90_up)-1)
for i in range(len(phi_13TeV_90_up)-1):
    avg_13TeV_90_up[i]=(phi_13TeV_90_up[i]+phi_13TeV_90_up[i+1])/2.
end = 2*phi_13TeV_90_up[-1]-avg_13TeV_90_up[-1]

avg_13TeV_90_up = np.append(avg_13TeV_90_up, end)

#적분값을 구해봅시다.
norm_90_up = 0.
for i in range(len(avg_13TeV_90_up)-1):
    if(-0.001<avg_13TeV_90_up[i]<np.pi):
        norm_90_up += data_13TeV_90_up[i+1]*(avg_13TeV_90_up[i+1]-avg_13TeV_90_up[i])

popt_13TeV_90_up, pcov_13TeV_90_up = curve_fit(wrap_templ_130(norm_90_up), phi_13TeV_90_up, data_13TeV_90_up, bounds=(0,100))

plt.scatter(phi_13TeV_90_up, data_13TeV_90_up, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_90_up))
G_periph0=np.zeros(len(phi_13TeV_90_up))
ridge_periph0=np.zeros(len(phi_13TeV_90_up))

wr.writerow([90,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_90_up)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_90_up[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_90_up[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_90_up[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_90_up, templ_130(phi_13TeV_90_up, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_90_up, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_90_up, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_90_up, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_90_up, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (90 \,\, \leq N^{rec}_{ch})$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_90~.png')

fig.clear()


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

# print(phi_13TeV_130_up)
# print(data_13TeV_130_up)
popt_13TeV_130_up, pcov_13TeV_130_up = curve_fit(wrap_templ_130(norm_130_up), phi_13TeV_130_up, data_13TeV_130_up, bounds=(0,100))

plt.scatter(phi_13TeV_130_up, data_13TeV_130_up, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_130_up))
G_periph0=np.zeros(len(phi_13TeV_130_up))
ridge_periph0=np.zeros(len(phi_13TeV_130_up))

wr.writerow([130,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_130_up)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_130_up[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_130_up, templ_130(phi_13TeV_130_up, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_130_up, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_130_up, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_130_up, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_130_up, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (130 \leq N^{rec}_{ch})$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_130~.png')

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

popt_13TeV_120_130, pcov_13TeV_120_130 = curve_fit(wrap_templ_130(norm_120_130), phi_13TeV_120_130, data_13TeV_120_130, bounds=(0,100))

plt.scatter(phi_13TeV_120_130, data_13TeV_120_130, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_120_130))
G_periph0=np.zeros(len(phi_13TeV_120_130))
ridge_periph0=np.zeros(len(phi_13TeV_120_130))

wr.writerow([125,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_120_130)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_120_130[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_120_130, templ_130(phi_13TeV_120_130, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_120_130, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_120_130, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_120_130, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_120_130, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (120 \leq N^{rec}_{ch} < 130)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_120~130.png')

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

popt_13TeV_110_120, pcov_13TeV_110_120 = curve_fit(wrap_templ_130(norm_110_120), phi_13TeV_110_120, data_13TeV_110_120, bounds=(0,100))

plt.scatter(phi_13TeV_110_120, data_13TeV_110_120, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_110_120))
G_periph0=np.zeros(len(phi_13TeV_110_120))
ridge_periph0=np.zeros(len(phi_13TeV_110_120))

wr.writerow([115,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_110_120)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_110_120[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_110_120, templ_130(phi_13TeV_110_120, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_110_120, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_110_120, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_110_120, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_110_120, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (110 \leq N^{rec}_{ch}<120)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_110~120.png')

fig.clear()


phi_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[0])
data_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[1])


norm = 0.

avg_13TeV_100_110=np.zeros(len(phi_13TeV_100_110)-1)
for i in range(len(phi_13TeV_100_110)-1):
    avg_13TeV_100_110[i]=(phi_13TeV_100_110[i]+phi_13TeV_100_110[i+1])/2.
end = 2*phi_13TeV_100_110[-1]-avg_13TeV_100_110[-1]

avg_13TeV_100_110 = np.append(avg_13TeV_100_110, end)

#적분값을 구해봅시다.
norm_100_110 = 0.
for i in range(len(avg_13TeV_100_110)-1):
    if(-0.001<avg_13TeV_100_110[i]<np.pi):
        norm_100_110 += data_13TeV_100_110[i+1]*(avg_13TeV_100_110[i+1]-avg_13TeV_100_110[i])

popt_13TeV_100_110, pcov_13TeV_100_110 = curve_fit(wrap_templ_130(norm_100_110), phi_13TeV_100_110, data_13TeV_100_110, bounds=(0,100))

plt.scatter(phi_13TeV_100_110, data_13TeV_100_110, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_100_110))
G_periph0=np.zeros(len(phi_13TeV_100_110))
ridge_periph0=np.zeros(len(phi_13TeV_100_110))

wr.writerow([105,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_100_110)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_100_110[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_100_110[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_100_110[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_100_110, templ_130(phi_13TeV_100_110, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_100_110, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_100_110, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_100_110, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_100_110, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (100 \leq N^{rec}_{ch}<110)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_100~110.png')

fig.clear()

phi_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[0])
data_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[1])

norm = 0.

avg_13TeV_90_100=np.zeros(len(phi_13TeV_90_100)-1)
for i in range(len(phi_13TeV_90_100)-1):
    avg_13TeV_90_100[i]=(phi_13TeV_90_100[i]+phi_13TeV_90_100[i+1])/2.
end = 2*phi_13TeV_90_100[-1]-avg_13TeV_90_100[-1]

avg_13TeV_90_100 = np.append(avg_13TeV_90_100, end)

#적분값을 구해봅시다.
norm_90_100 = 0.
for i in range(len(avg_13TeV_90_100)-1):
    if(-0.001<avg_13TeV_90_100[i]<np.pi):
        norm_90_100 += data_13TeV_90_100[i+1]*(avg_13TeV_90_100[i+1]-avg_13TeV_90_100[i])

popt_13TeV_90_100, pcov_13TeV_90_100 = curve_fit(wrap_templ_130(norm_90_100), phi_13TeV_90_100, data_13TeV_90_100, bounds=(0,100))

plt.scatter(phi_13TeV_90_100, data_13TeV_90_100, color = 'black', s=1000, marker='o',label=fr'$Y(\Delta\phi)$')

G_periph=np.zeros(len(phi_13TeV_90_100))
G_periph0=np.zeros(len(phi_13TeV_90_100))
ridge_periph0=np.zeros(len(phi_13TeV_90_100))

wr.writerow([95,Ffit[-1],Gfit[-1],vfit[-1],'13TeV'])

for i in range(len(phi_13TeV_90_100)):
    G_periph[i] = templ_130_G_periph(phi_13TeV_90_100[i], Gfit[-1], Ffit[-1], vfit[-1], i)
    G_periph0[i] = templ_130_G_periph0(phi_13TeV_90_100[i], Gfit[-1], Ffit[-1], vfit[-1])
    ridge_periph0[i] = templ_130_ridge_periph0(phi_13TeV_90_100[i], Gfit[-1], Ffit[-1], vfit[-1])

plt.plot(phi_13TeV_90_100, templ_130(phi_13TeV_90_100, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')

# plt.plot(phi_13TeV_90_100, G_periph, drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$G+FY^periph(\Delta\phi)$')
plt.scatter(phi_13TeV_90_100, G_periph, edgecolor = 'black', facecolors='none', s=1000, label=r'$G+FY^{periph} (\Delta\phi)$')
plt.plot(phi_13TeV_90_100, G_periph0, color = 'orange', linewidth=5, linestyle = '--',label=r'$G+FY^{periph}(0)$')
plt.plot(phi_13TeV_90_100, ridge_periph0, color = 'blue', linewidth=5, linestyle = '--',label=r'$Y^{ridge}+FY^{periph}(0)$')

plt.xlim(-0.96,0.96)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS \,\, pp(13TeV) \,\, (90 \,\, \leq N^{rec}_{ch}<100)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_90~100.png')

fig.clear()


ridge_parameters.close()        #다시 불러들이기 위해 파일을 저장하자. 
                                #이 그래프 이전에는 데이터 피팅, 이 그래프 이후에는 fitting한 F, G, v를 이용한 그래프를 따로 그려야 함.
v22_result_multi=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[0],skiprows=1)
v22_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[3],skiprows=1)
v22_data_multi=np.loadtxt('./atlasgraphs/v22_13TeV.csv',delimiter=',',usecols=[0])
v22_data=np.loadtxt('./atlasgraphs/v22_13TeV.csv',delimiter=',',usecols=[1])

plt.scatter(v22_result_multi, v22_result, color = 'black', s=1000, label=r'$result$')
plt.scatter(v22_data_multi, v22_data, facecolor = 'none', edgecolor = 'red', linewidths=5, s=1000, label=r'$data\,thief$')
plt.legend(fontsize=45,framealpha=False)

plt.scatter(v22_data_multi, v22_data, color = 'red', marker='+', linewidths=5, s=1000, label=r'$data\,thief$')

plt.ylabel(r'$v_{2,2}$', size = 70)
plt.xlabel(r'$N_{ch}^{rec}$', size=70)

plt.minorticks_on()

plt.ylim(0,0.0073)
plt.xlim(15,140)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/v22_13TeV.png')

fig.clear()


G_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[2],skiprows=1)
atlas_phi = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[3],skiprows=1)
#계산에 사용한 atlas delta phi는 모두 동일하므로 통일해도 관계없다.

atlas_phi_95 = np.loadtxt('phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_95 = np.loadtxt('phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_105 = np.loadtxt('phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_105 = np.loadtxt('phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_115 = np.loadtxt('phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_115 = np.loadtxt('phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_125 = np.loadtxt('phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_125 = np.loadtxt('phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_135 = np.loadtxt('phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_135 = np.loadtxt('phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[3],skiprows=1)



def atlas_ridge(phi, G, v):          #atlas는 0.5<pt<5 이므로 지금까지 그려온 그래프와 맞지 않음.
        return G*(1+2*v*np.cos(2*phi))

Deltaphi = np.arange(-0.96,0.96,0.04)
Deltaphi = np.append(Deltaphi,0.96)

# for i in range(4):    Y^ridge, for Experiment data
Multi_90_up = atlas_ridge(Deltaphi, G_result[0], v22_result[0])
Multi_130_up = atlas_ridge(Deltaphi, G_result[1], v22_result[1])
Multi_120_130 = atlas_ridge(Deltaphi, G_result[2], v22_result[2])
Multi_110_120 = atlas_ridge(Deltaphi, G_result[3], v22_result[3])
Multi_100_110 = atlas_ridge(Deltaphi, G_result[4], v22_result[4])
Multi_90_100 = atlas_ridge(Deltaphi, G_result[5], v22_result[5])

# data-Czyam의 적분값 구하기.
norm_templ_135 = 0.
norm_templ_125 = 0.
norm_templ_115 = 0.
norm_templ_105 = 0.
norm_templ_95 = 0.

for i in range(len(data_13TeV_130_up)):
    norm_templ_135 += (data_13TeV_130_up[i]-min(data_13TeV_130_up))*(avg_13TeV_130_up[i]-avg_13TeV_130_up[i-1])
    norm_templ_125 += (data_13TeV_120_130[i]-min(data_13TeV_120_130))*(avg_13TeV_120_130[i]-avg_13TeV_120_130[i-1])
    norm_templ_115 += (data_13TeV_110_120[i]-min(data_13TeV_110_120))*(avg_13TeV_110_120[i]-avg_13TeV_110_120[i-1])
    norm_templ_105 += (data_13TeV_100_110[i]-min(data_13TeV_100_110))*(avg_13TeV_100_110[i]-avg_13TeV_100_110[i-1])
    norm_templ_95 += (data_13TeV_90_100[i]-min(data_13TeV_90_100))*(avg_13TeV_90_100[i]-avg_13TeV_90_100[i-1])
    # print(data_13TeV_130_up[i]-min(data_13TeV_130_up), avg_13TeV_130_up[i], avg_13TeV_130_up[i-1])

# print(norm_templ_135)
theory_norm = np.loadtxt('phiCorr_norm.csv', delimiter=',', usecols=[1])

# print(theory_norm)

#그냥 뺄 경우 theory-periph에 나타나는 음수를 주의하자.
theory_norm[0] -= 2*0.96*min(atlas_result_135)
theory_norm[1] -= 2*0.96*min(atlas_result_125)
theory_norm[2] -= 2*0.96*min(atlas_result_115)
theory_norm[3] -= 2*0.96*min(atlas_result_105)
theory_norm[4] -= 2*0.96*min(atlas_result_95)

# print(theory_norm)

ratio_95 = norm_templ_95/theory_norm[4]
ratio_105 = norm_templ_105/theory_norm[3]
ratio_115 = norm_templ_115/theory_norm[2]
ratio_125 = norm_templ_125/theory_norm[1]
ratio_135 = norm_templ_135/theory_norm[0]

print(ratio_135)

theory_norm_95 = (atlas_result_95-min(atlas_result_95))*ratio_95
theory_norm_105 = (atlas_result_105-min(atlas_result_105))*ratio_105
theory_norm_115 = (atlas_result_115-min(atlas_result_115))*ratio_115
theory_norm_125 = (atlas_result_125-min(atlas_result_125))*ratio_125
theory_norm_135 = (atlas_result_135-min(atlas_result_135))*ratio_135

fig.set_size_inches(40, 20, forward=True)

#단순 theory 계산값 -Czyam
# plt.plot(atlas_phi_95, atlas_result_95-min(atlas_result_95), color = 'blue', linewidth=5, linestyle = '-',label=r'$result$')
plt.plot(atlas_phi_135, theory_norm_135, color = 'black', linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi_125, theory_norm_125, color = 'green', linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi_115, theory_norm_115, color = 'orange', linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi_105, theory_norm_105, color = 'blue', linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
# plt.plot(atlas_phi_95, theory_norm_95, color = 'grey', linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')
# plt.plot(atlas_phi_135, atlas_result_135-min(atlas_result_135), color = 'black', linewidth=5, linestyle = ':',label=r'$130 \leq N_{ch}^{rec}, ratio x$')

plt.legend(fontsize=45,framealpha=False,loc='upper right')

#atlas의 Yridge를 계산한 값 -Czyam
# plt.plot(Deltaphi, Multi_90_up-min(Multi_90_up), color = 'red', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_130_up-min(Multi_130_up), color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_120_130-min(Multi_120_130), color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
# plt.plot(Deltaphi, Multi_110_120-min(Multi_110_120), color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
# plt.plot(Deltaphi, Multi_100_110-min(Multi_100_110), color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
# plt.plot(Deltaphi, Multi_90_100-min(Multi_90_100), color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')

# print(Deltaphi)


#단순히 data-Czyam
# plt.plot(phi_13TeV_90_up, data_13TeV_90_up-min(data_13TeV_90_up), marker = 'o', markersize = 15, color = 'red', linewidth=5, linestyle = ':',label=r'$90 \,\, \leq N_{ch}^{rec}$')
plt.plot(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), marker = 'o', markersize = 15, color = 'black', linewidth=5, linestyle = ':',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), marker = 'o', markersize = 15, color = 'green', linewidth=5, linestyle = ':',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), marker = 'o', markersize = 15, color = 'orange', linewidth=5, linestyle = ':',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), marker = 'o', markersize = 15, color = 'blue', linewidth=5, linestyle = ':',label=r'$100 \leq N_{ch}^{rec}<110$')
# plt.plot(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), marker = 'o', markersize = 15, color = 'grey', linewidth=5, linestyle = ':',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')



plt.ylabel(r'$Y^{data}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.,1.)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.13,0.4), ncol=2)
# plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi-czyam.png')

fig.clear()


plt.plot(atlas_phi, atlas_result, color = 'blue', linewidth=5, linestyle = '-',label=r'$result$')
plt.plot(Deltaphi, Multi_130_up, color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_120_130, color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(Deltaphi, Multi_110_120, color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(Deltaphi, Multi_100_110, color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')

plt.ylabel(r'$Y^{ridge}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.,1.)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))
plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi_Yridge.png')

fig.clear()

# theory - Czyam    //    Yridge - Czyam 비교

F_result_multi=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[0],skiprows=1)
F_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[1],skiprows=1)

#Theory Result
# theory_90_up = atlas_result-F_result[0]*periph(atlas_phi, *popt_periph_130TeV)
theory_130_up = theory_norm_135-F_result[1]*periph(atlas_phi_135, *popt_periph_130TeV)
theory_120_130 = theory_norm_125-F_result[2]*periph(atlas_phi_125, *popt_periph_130TeV)
theory_110_120 = theory_norm_115-F_result[3]*periph(atlas_phi_115, *popt_periph_130TeV)
theory_100_110 = theory_norm_105-F_result[4]*periph(atlas_phi_105, *popt_periph_130TeV)
theory_90_100 = theory_norm_95-F_result[5]*periph(atlas_phi_95, *popt_periph_130TeV)

# print(F_result)

# plt.plot(atlas_phi, theory_90_up-min(theory_90_up), color = 'red', linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}$')
plt.plot(atlas_phi_135, theory_130_up-min(theory_130_up), color = 'black', linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi_125, theory_120_130-min(theory_120_130), color = 'green', linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi_115, theory_110_120-min(theory_110_120), color = 'orange', linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi_105, theory_100_110-min(theory_100_110), color = 'blue', linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(atlas_phi_95, theory_90_100-min(theory_90_100), color = 'grey', linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')

plt.legend(fontsize=45,framealpha=False,loc='upper right')

#experiment result
# plt.plot(Deltaphi, Multi_90_up-min(Multi_90_up), color = 'red', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_130_up-min(Multi_130_up), color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_120_130-min(Multi_120_130), color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(Deltaphi, Multi_110_120-min(Multi_110_120), color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(Deltaphi, Multi_100_110-min(Multi_100_110), color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(Deltaphi, Multi_90_100-min(Multi_90_100), color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')


plt.ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.,1.)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))



plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi_Yridge_czyam.png')

fig.clear()

#C_zyam 제외

theory_90_up = atlas_result-F_result[0]*periph(atlas_phi, *popt_periph_130TeV)
theory_130_up = atlas_result-F_result[1]*periph(atlas_phi, *popt_periph_130TeV)
theory_120_130 = atlas_result-F_result[2]*periph(atlas_phi, *popt_periph_130TeV)
theory_110_120 = atlas_result-F_result[3]*periph(atlas_phi, *popt_periph_130TeV)
theory_100_110 = atlas_result-F_result[4]*periph(atlas_phi, *popt_periph_130TeV)
theory_90_100 = atlas_result-F_result[5]*periph(atlas_phi, *popt_periph_130TeV)

plt.plot(atlas_phi, theory_90_up, color = 'red', linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_130_up, color = 'black', linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_120_130, color = 'green', linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi, theory_110_120, color = 'orange', linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi, theory_100_110, color = 'blue', linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(atlas_phi, theory_90_100, color = 'grey', linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')

plt.legend(fontsize=45,framealpha=False)

#experiment result
# plt.plot(Deltaphi, Multi_90_up color = 'red', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_130_up, color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_120_130, color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
# plt.plot(Deltaphi, Multi_110_120, color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
# plt.plot(Deltaphi, Multi_100_110, color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
# plt.plot(Deltaphi, Multi_90_100, color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')


plt.ylabel(r'$Y^{ridge}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.,1.)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))



plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_Ytheory.png')

fig.clear()