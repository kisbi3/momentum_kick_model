import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv


mpl.rcParams["text.usetex"] = True
fig = plt.figure()
ax = plt.axes()
fig.set_size_inches(25, 20, forward=True)

rapid_000_rapid=np.loadtxt('./rapidity/hdndy_mult000.txt',delimiter='\t',usecols=[0])
rapid_000_yield=np.loadtxt('./rapidity/hdndy_mult000.txt',delimiter='\t',usecols=[1])
rapid_000_ynorm=np.loadtxt('./rapidity/hdndy_mult000.txt',delimiter='\t',usecols=[2])
rapid_035_rapid=np.loadtxt('./rapidity/hdndy_mult035.txt',delimiter='\t',usecols=[0])
rapid_035_yield=np.loadtxt('./rapidity/hdndy_mult035.txt',delimiter='\t',usecols=[1])
rapid_035_ynorm=np.loadtxt('./rapidity/hdndy_mult035.txt',delimiter='\t',usecols=[2])
rapid_080_rapid=np.loadtxt('./rapidity/hdndy_mult080.txt',delimiter='\t',usecols=[0])
rapid_080_yield=np.loadtxt('./rapidity/hdndy_mult080.txt',delimiter='\t',usecols=[1])
rapid_080_ynorm=np.loadtxt('./rapidity/hdndy_mult080.txt',delimiter='\t',usecols=[2])
rapid_105_rapid=np.loadtxt('./rapidity/hdndy_mult105.txt',delimiter='\t',usecols=[0])
rapid_105_yield=np.loadtxt('./rapidity/hdndy_mult105.txt',delimiter='\t',usecols=[1])
rapid_105_ynorm=np.loadtxt('./rapidity/hdndy_mult105.txt',delimiter='\t',usecols=[2])

# print(rapid_000_rapid)
# print(rapid_000_yield)

plt.plot(rapid_000_rapid, rapid_000_yield-min(rapid_000_yield), color = 'red', linewidth=5, linestyle = '-', label=fr'$total$')
plt.plot(rapid_035_rapid, rapid_035_yield-min(rapid_035_yield), color = 'green', linewidth=5, linestyle = '-', label=fr'$N \geq 35$')
plt.plot(rapid_080_rapid, rapid_080_yield-min(rapid_080_yield), color = 'blue', linewidth=5, linestyle = '-', label=fr'$N \geq 85$')
plt.plot(rapid_105_rapid, rapid_105_yield-min(rapid_105_yield), color = 'black', linewidth=5, linestyle = '-', label=fr'$N \geq 105$')

plt.ylabel(r'$Y$', size = 70)
plt.xlabel(r'$y(rapidity)$', size=70)

# plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()
plt.legend(fontsize=45,framealpha=False,loc='upper right')
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/rapidity.png')

fig.clear()

plt.plot(rapid_000_rapid, rapid_000_ynorm-min(rapid_000_ynorm), color = 'red', linewidth=5, linestyle = '-', label=fr'$total$')
plt.plot(rapid_035_rapid, rapid_035_ynorm-min(rapid_035_ynorm), color = 'green', linewidth=5, linestyle = '-', label=fr'$N \geq 35$')
plt.plot(rapid_080_rapid, rapid_080_ynorm-min(rapid_080_ynorm), color = 'blue', linewidth=5, linestyle = '-', label=fr'$N \geq 85$')
plt.plot(rapid_105_rapid, rapid_105_ynorm-min(rapid_105_ynorm), color = 'black', linewidth=5, linestyle = '-', label=fr'$N \geq 105$')

plt.ylabel(r'$Y$', size = 70)
plt.xlabel(r'$y(rapidity)$', size=70)

# plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()
plt.legend(fontsize=45,framealpha=False,loc='upper right')
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/rapidity_norm.png')

fig.clear()


pt_000_ptdis=np.loadtxt('./pt/hdndpt_mult000.txt',delimiter='\t',usecols=[0])
pt_000_yield=np.loadtxt('./pt/hdndpt_mult000.txt',delimiter='\t',usecols=[1])
pt_000_ynorm=np.loadtxt('./pt/hdndpt_mult000.txt',delimiter='\t',usecols=[2])
pt_035_ptdis=np.loadtxt('./pt/hdndpt_mult035.txt',delimiter='\t',usecols=[0])
pt_035_yield=np.loadtxt('./pt/hdndpt_mult035.txt',delimiter='\t',usecols=[1])
pt_035_ynorm=np.loadtxt('./pt/hdndpt_mult035.txt',delimiter='\t',usecols=[2])
pt_080_ptdis=np.loadtxt('./pt/hdndpt_mult080.txt',delimiter='\t',usecols=[0])
pt_080_yield=np.loadtxt('./pt/hdndpt_mult080.txt',delimiter='\t',usecols=[1])
pt_080_ynorm=np.loadtxt('./pt/hdndpt_mult080.txt',delimiter='\t',usecols=[2])
pt_105_ptdis=np.loadtxt('./pt/hdndpt_mult105.txt',delimiter='\t',usecols=[0])
pt_105_yield=np.loadtxt('./pt/hdndpt_mult105.txt',delimiter='\t',usecols=[1])
pt_105_ynorm=np.loadtxt('./pt/hdndpt_mult105.txt',delimiter='\t',usecols=[2])

# print(rapid_000_rapid)
# print(rapid_000_yield)

plt.plot(pt_000_ptdis, pt_000_yield-min(pt_000_yield), color = 'red', linewidth=5, linestyle = '-', label=fr'$total$')
plt.plot(pt_035_ptdis, pt_035_yield-min(pt_035_yield), color = 'green', linewidth=5, linestyle = '-', label=fr'$N \geq 35$')
plt.plot(pt_080_ptdis, pt_080_yield-min(pt_080_yield), color = 'blue', linewidth=5, linestyle = '-', label=fr'$N \geq 85$')
plt.plot(pt_105_ptdis, pt_105_yield-min(pt_105_yield), color = 'black', linewidth=5, linestyle = '-', label=fr'$N \geq 105$')

plt.ylabel(r'$Y$', size = 70)
plt.xlabel(r'$p_T$', size=70)

plt.xlim(0,3)

# plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()
plt.legend(fontsize=45,framealpha=False,loc='upper right')
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/pt.png')

fig.clear()


plt.plot(pt_000_ptdis, pt_000_ynorm-min(pt_000_ynorm), color = 'red', linewidth=5, linestyle = '-', label=fr'$total$')
plt.plot(pt_035_ptdis, pt_035_ynorm-min(pt_035_ynorm), color = 'green', linewidth=5, linestyle = '-', label=fr'$N \geq 35$')
plt.plot(pt_080_ptdis, pt_080_ynorm-min(pt_080_ynorm), color = 'blue', linewidth=5, linestyle = '-', label=fr'$N \geq 85$')
plt.plot(pt_105_ptdis, pt_105_ynorm-min(pt_105_ynorm), color = 'black', linewidth=5, linestyle = '-', label=fr'$N \geq 105$')

plt.ylabel(r'$Y$', size = 70)
plt.xlabel(r'$p_T$', size=70)

plt.xlim(0,3)

# plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()
plt.legend(fontsize=45,framealpha=False,loc='upper right')
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/pt_norm.png')

fig.clear()