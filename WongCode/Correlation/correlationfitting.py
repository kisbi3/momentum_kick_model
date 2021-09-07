import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

deltaphi=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphi=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
data_error1=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
data_error2=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('Correlation_a03.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi1=np.loadtxt('Correlation_a03.csv',delimiter=',',usecols=[1],skiprows=1)
resultdNdphi2=np.loadtxt('Correlation_a05.csv',delimiter=',',usecols=[1],skiprows=1)
resultdNdphi3=np.loadtxt('Correlation_a07.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin1=min(resultdNdphi1)
dNdphimin2=min(resultdNdphi2)
dNdphimin3=min(resultdNdphi3)
resultdNdphi_removeE=np.loadtxt('Correlation_removeE.csv',delimiter=',',usecols=[1],skiprows=1)

# print(dNdphimin)


plt.errorbar(deltaphi,dNdphi, yerr=(data_error1,abs(data_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \,$CMS',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
plt.plot(deltaphi,dNdphi, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, resultdNdphi1-dNdphimin1, color = 'red', linewidth=7, linestyle = '-',label=r'$a=0.3$')
plt.plot(resultphi, resultdNdphi2-dNdphimin2, color = 'green', linewidth=7, linestyle = '-',label=r'$a=0.5$')
plt.plot(resultphi, resultdNdphi3-dNdphimin3, color = 'blue', linewidth=7, linestyle = '-',label=r'$a=0.7$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Correlation_fitting.png')