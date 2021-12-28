import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import pandas as pd


# f = open('G_coefficient.csv','w',newline='')
# wr=csv.writer(f)

# df = pd.read_csv('G_coefficient.csv', delimiter=',')
dataintegral=[]     #data 적분 결과 저장 array

mpl.rcParams["text.usetex"] = True

def periph(phi, a, b):
    return a*np.cos(phi)+b


def templ_276(phi, dataintegral, F, v):
    G=(dataintegral-F*norm_276)/np.pi
    return F*periph_276_data+G*(1+2*v*np.cos(2*phi))


def templ_130_graph(phi, G, F, v):
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
        # wr.writerow([G])
        # G = np.append(G, G_dist)\
        print(G_dist, F, v)
        Gfit.append(G_dist)
        Ffit.append(F)
        vfit.append(v)
        return templ_130_graph(phi, G_dist, F, v)
    
    # print(Gfit)
    return tempfunc
        # return F*periph_130_data+G(1+2*v*np.cos(2*phi))


fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

periph_276_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[0],skiprows=1)
periph_276_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[1],skiprows=1)
periph_130_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_130_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)

# plt.plot(periph_276_phi, periph_276_data, color = 'red', linewidth=7, linestyle = '-',label=fr'$peripheral, 2.76TeV$')
# plt.plot(periph_130_phi, periph_130_data, color = 'blue', linewidth=7, linestyle = '-',label=fr'$peripheral, 13TeV$')

plt.scatter(periph_276_phi, periph_276_data, color = 'red', s=100,label=fr'$peripheral, 2.76TeV$')
plt.scatter(periph_130_phi, periph_130_data, color = 'blue', s=100,label=fr'$peripheral, 13TeV$')

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)

# x = np.linspace(-1.5, 5, 35)

popt_276, pcov_276 = curve_fit(periph, periph_276_phi, periph_276_data)
popt_130, pcov_130 = curve_fit(periph, periph_130_phi, periph_130_data)

# plt.plot(periph_276_phi, periph(periph_276_phi, *popt_276), color = 'red', linewidth=7, linestyle = '-',label=fr'$peripheral\,\, fit, 2.76TeV$')
# plt.plot(periph_130_phi, periph(periph_130_phi, *popt_130), color = 'blue', linewidth=7, linestyle = '-',label=fr'$peripheral\,\, fit, 13TeV$')

# print(popt_130, popt_276)


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

# print(norm_130, norm_276)
# print(periph_130_phi_avg)


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

print(norm_130_up)
dataintegral.append(norm_130_up)
        # print(data_13TeV_130_up[i+1], avg_13TeV_130_up[i+1], avg_13TeV_130_up[i])
    # if(-0.001<periph_276_phi_avg[i]<np.pi):
    #     norm_276 += periph_276_data[i+1]*(periph_276_phi_avg[i+1]-periph_276_phi_avg[i])

# wr.writerow(['130~'])


popt_13TeV_130_up, pcov_13TeV_130_up = curve_fit(wrap_templ_130(norm_130_up), phi_13TeV_130_up, data_13TeV_130_up)
# popt_13TeV_130_up, pcov_13TeV_130_up = curve_fit(templ_130, phi_13TeV_130_up, data_13TeV_130_up)
# print(G[-1])

print(popt_13TeV_130_up)
# print(pcov_13TeV_130_up)

# print(Gfit[-1])
# print(Ffit[-1])
# print(vfit[-1])


# plt.scatter(phi_13TeV_130_up, templ_130(phi_13TeV_130_up, *popt_13TeV_130_up), color = 'red', s=5000,marker='_',label=fr'$templ, 13TeV$')
# plt.step(phi_13TeV_130_up, templ_130(phi_13TeV_130_up, *popt_13TeV_130_up), drawstyle='steps-mid', color = 'red', label=fr'$templ, 13TeV$')
plt.scatter(phi_13TeV_130_up, data_13TeV_130_up, color = 'blue', s=1000, marker='o',label=fr'$data, 13TeV$')
# phi_13TeV_130_up_plus = np.append(phi_13TeV_130_up, -1.55)
# plt.plot(phi_13TeV_130_up, templ_130(phi_13TeV_130_up, *popt_13TeV_130_up), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$temple\,\, fit, 13TeV$')
plt.plot(phi_13TeV_130_up, templ_130_graph(phi_13TeV_130_up, Gfit[-1], Ffit[-1], vfit[-1]), drawstyle='steps-mid', color = 'red', linewidth=5, linestyle = '-',label=fr'$temple\,\, fit, 13TeV$')

plt.xlim(-1.55,4.7)

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS,\quad Y^{templ} \, (130 \leq N^{rec}_{ch})$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.legend(fontsize=45,framealpha=False)


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('13TeV_130~.png')

fig.clear()



