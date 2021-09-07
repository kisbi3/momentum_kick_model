import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import csv


# mpl.rcParams["text.usetex"] = True

fig = plt.figure(figsize = (25,18))
ax = plt.axes()

fig.set_size_inches(23.384, 18, forward=True)

# pt=np.loadtxt('Test.csv',delimiter=',',usecols=[0],skiprows=1)
# y0=np.loadtxt('Test.csv',delimiter=',',usecols=[1],skiprows=1)
# y1=np.loadtxt('Test.csv',delimiter=',',usecols=[2],skiprows=1)
# y2=np.loadtxt('Test.csv',delimiter=',',usecols=[3],skiprows=1)
# y3=np.loadtxt('Test.csv',delimiter=',',usecols=[4],skiprows=1)
# y4=np.loadtxt('Test.csv',delimiter=',',usecols=[5],skiprows=1)

# # pt=np.loadtxt('Test_3.csv',delimiter=',',usecols=[0],skiprows=1)
# # y30=np.loadtxt('Test_3.csv',delimiter=',',usecols=[1],skiprows=1)
# # y31=np.loadtxt('Test_3.csv',delimiter=',',usecols=[2],skiprows=1)
# # y32=np.loadtxt('Test_3.csv',delimiter=',',usecols=[3],skiprows=1)
# # y33=np.loadtxt('Test_3.csv',delimiter=',',usecols=[4],skiprows=1)
# # y34=np.loadtxt('Test_3.csv',delimiter=',',usecols=[5],skiprows=1)

# Refpt0=np.loadtxt('data/pT_fixed_y_0.csv',delimiter=',',usecols=[0])
# Refy0=np.loadtxt('data/pT_fixed_y_0.csv',delimiter=',',usecols=[1])
# Refpt1=np.loadtxt('data/pT_fixed_y_1.csv',delimiter=',',usecols=[0])
# Refy1=np.loadtxt('data/pT_fixed_y_1.csv',delimiter=',',usecols=[1])
# Refpt2=np.loadtxt('data/pT_fixed_y_2.csv',delimiter=',',usecols=[0])
# Refy2=np.loadtxt('data/pT_fixed_y_2.csv',delimiter=',',usecols=[1])
# Refpt3=np.loadtxt('data/pT_fixed_y_3.csv',delimiter=',',usecols=[0])
# Refy3=np.loadtxt('data/pT_fixed_y_3.csv',delimiter=',',usecols=[1])
# Refpt4=np.loadtxt('data/pT_fixed_y_4.csv',delimiter=',',usecols=[0])
# Refy4=np.loadtxt('data/pT_fixed_y_4.csv',delimiter=',',usecols=[1])

# plt.xlabel(r'$p_{t}\quad(GeV/c)$',size=50)
# plt.ylabel(r'$dF/dyp_{t}dp_{t}\quad(GeV^{-2})$',size=50)
# plt.plot(pt,y0, color="black",linewidth=3,label=r'$y=0$', linestyle = '-')
# plt.plot(pt,y1, color="red",linewidth=3,label=r'$y=1$', linestyle = '-')
# plt.plot(pt,y2, color="blue",linewidth=3,label=r'$y=2$', linestyle = '-')
# plt.plot(pt,y3, color="magenta",linewidth=3,label=r'$y=3$', linestyle = '-')
# plt.plot(pt,y4, color="green",linewidth=3,label=r'$y=4$', linestyle = '-')

# # plt.plot(pt,y30, color="black",linewidth=3,label=r'$y=0$', linestyle = ':')
# # plt.plot(pt,y31, color="red",linewidth=3,label=r'$y=1$', linestyle = ':')
# # plt.plot(pt,y32, color="blue",linewidth=3,label=r'$y=2$', linestyle = ':')
# # plt.plot(pt,y33, color="magenta",linewidth=3,label=r'$y=3$', linestyle = ':')
# # plt.plot(pt,y34, color="green",linewidth=3,label=r'$y=4$', linestyle = ':')


# # plt.scatter(Refpt0,Refy0, s = 200, c = 'black')
# # plt.scatter(Refpt1,Refy1, s = 200, c = 'red')
# # plt.scatter(Refpt2,Refy2, s = 200, c = 'blue')
# # plt.scatter(Refpt3,Refy3, s = 200, c = 'magenta')
# # plt.scatter(Refpt4,Refy4, s = 200, c = 'green')

# plt.minorticks_on()
# plt.yscale('log')
# ax.axis([0,4,0.00001,2.5]) #2.5
# # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


# plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False)
# plt.tight_layout()

# fig.savefig('Fix_y_pivs3pi.png')
# # plt.show()

# fig.clear()

# yi=np.loadtxt('Test2.csv',delimiter=',',usecols=[0],skiprows=1)
# pt0=np.loadtxt('Test2.csv',delimiter=',',usecols=[1],skiprows=1)
# pt1=np.loadtxt('Test2.csv',delimiter=',',usecols=[2],skiprows=1)
# pt2=np.loadtxt('Test2.csv',delimiter=',',usecols=[3],skiprows=1)
# pt3=np.loadtxt('Test2.csv',delimiter=',',usecols=[4],skiprows=1)
# pt4=np.loadtxt('Test2.csv',delimiter=',',usecols=[5],skiprows=1)

# # pt30=np.loadtxt('Test2_3.csv',delimiter=',',usecols=[1],skiprows=1)
# # pt31=np.loadtxt('Test2_3.csv',delimiter=',',usecols=[2],skiprows=1)
# # pt32=np.loadtxt('Test2_3.csv',delimiter=',',usecols=[3],skiprows=1)
# # pt33=np.loadtxt('Test2_3.csv',delimiter=',',usecols=[4],skiprows=1)
# # pt34=np.loadtxt('Test2_3.csv',delimiter=',',usecols=[5],skiprows=1)

# Refy00=np.loadtxt('data/y_fixed_pT_0.csv',delimiter=',',usecols=[0])
# Refpt00=np.loadtxt('data/y_fixed_pT_0.csv',delimiter=',',usecols=[1])
# Refy01=np.loadtxt('data/y_fixed_pT_1.csv',delimiter=',',usecols=[0])
# Refpt01=np.loadtxt('data/y_fixed_pT_1.csv',delimiter=',',usecols=[1])
# Refy02=np.loadtxt('data/y_fixed_pT_2.csv',delimiter=',',usecols=[0])
# Refpt02=np.loadtxt('data/y_fixed_pT_2.csv',delimiter=',',usecols=[1])
# Refy03=np.loadtxt('data/y_fixed_pT_3.csv',delimiter=',',usecols=[0])
# Refpt03=np.loadtxt('data/y_fixed_pT_3.csv',delimiter=',',usecols=[1])
# Refy04=np.loadtxt('data/y_fixed_pT_4.csv',delimiter=',',usecols=[0])
# Refpt04=np.loadtxt('data/y_fixed_pT_4.csv',delimiter=',',usecols=[1])

# plt.xlabel(r'$y$',size=50)
# plt.ylabel(r'$dF/dyp_{t}dp_{t}\quad(GeV^{-2})$',size=50)
# plt.plot(yi,pt0, color="black",linewidth=3,label=r'$pt=0.1$', linestyle = '-')
# plt.plot(yi,pt1, color="red",linewidth=3,label=r'$pt=1.$', linestyle = '-')
# plt.plot(yi,pt2, color="blue",linewidth=3,label=r'$pt=2$', linestyle = '-')
# plt.plot(yi,pt3, color="magenta",linewidth=3,label=r'$pt=3$', linestyle = '-')
# plt.plot(yi,pt4, color="green",linewidth=3,label=r'$pt=4$', linestyle = '-')



# # plt.plot(yi,pt30, color="black",linewidth=3,label=r'$pt=0$', linestyle = ':')
# # plt.plot(yi,pt31, color="red",linewidth=3,label=r'$pt=1$', linestyle = ':')
# # plt.plot(yi,pt32, color="blue",linewidth=3,label=r'$pt=2$', linestyle = ':')
# # plt.plot(yi,pt33, color="magenta",linewidth=3,label=r'$pt=3$', linestyle = ':')
# # plt.plot(yi,pt34, color="green",linewidth=3,label=r'$pt=4$', linestyle = ':')

# # plt.scatter(Refy00,Refpt00, s = 200, c = 'black')
# # plt.scatter(Refy01,Refpt01, s = 200, c = 'red')
# # plt.scatter(Refy02,Refpt02, s = 200, c = 'blue')
# # plt.scatter(Refy03,Refpt03, s = 200, c = 'magenta')
# # plt.scatter(Refy04,Refpt04, s = 200, c = 'green')


# plt.minorticks_on()
# plt.yscale('log')
# # ax.axis([0,4,0.00004,10])
# plt.xlim(-6,6)
# plt.ylim(0.00001,2.5) #2.5
# # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# # ax.spines['left'].set_position('center')
# # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


# plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.tight_layout()

# fig.savefig('Fix_pt_pivs3pi.png')

# fig.clear()

yi=np.loadtxt('dNdy.csv',delimiter=',',usecols=[0],skiprows=1)
dNdy=np.loadtxt('dNdy.csv',delimiter=',',usecols=[1],skiprows=1)

# plt.scatter(yi,dNdy, s = 200, c = 'black')
plt.plot(yi,dNdy, color="black",linewidth=5, linestyle = '-')


plt.minorticks_on()
plt.yscale('log')
plt.ylim(0.000001,10)
plt.xlim(0.0,5.5)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# ax.spines['left'].set_position('center')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('dFdy.png')