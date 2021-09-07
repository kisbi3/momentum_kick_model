import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from multiprocessing import Process
from matplotlib import cm
import os


mpl.rcParams["text.usetex"] = True

fig = plt.figure()
fig.set_size_inches(15, 10, forward=True)

detaf=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[0],skiprows=1)
dphif=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)



def functions(numbers):
    proc = os.getpid()
    if numbers == 1:

        pt015=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[2],skiprows=1)
        ax3d = plt.axes(projection="3d")
        ax3d.plot_trisurf(dphif,detaf,pt015,cmap='jet', linewidths=0.5)

        plt.title(r'$p_{T}=0.15,Ridge$', size = 40)
        ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
        ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
        ax3d.set_zlabel(r"$dN/{d\Delta \eta d\Delta\phi p_{t}dp_{t}}$", size=25, labelpad = 25)
        # axes.set_zlabel("$\Delta\eta$")
        plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
        plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

        plt.tight_layout()
        # plt.show()

        fig.savefig('pt=015_Ridge_3dplot.png')
        fig.clear()

    if numbers == 2:

        ax3d = plt.axes(projection="3d")
        pt1=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[3],skiprows=1)
        ax3d.plot_trisurf(dphif,detaf,pt1,cmap='jet', linewidths=0.5)

        plt.title(r'$p_{T}=1,Ridge$', size = 40)
        ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
        ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
        # axes.set_zlabel("$\Delta\eta$")
        plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
        plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

        plt.tight_layout()
        # plt.show()

        fig.savefig('pt=1_Ridge_3dplot.png')

        fig.clear()

    if numbers == 3:
        ax3d = plt.axes(projection="3d")
        pt2=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[4],skiprows=1)
        ax3d.plot_trisurf(dphif,detaf,pt2,cmap='jet', linewidths=0.5)

        plt.title(r'$p_{T}=2,Ridge$', size = 40)
        ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
        ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
        # axes.set_zlabel("$\Delta\eta$")
        plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
        plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

        plt.tight_layout()
        # plt.show()

        fig.savefig('pt=2_Ridge_3dplot.png')

        fig.clear()
    if numbers == 4:
        ax3d = plt.axes(projection="3d")
        pt3=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[5],skiprows=1)
        ax3d.plot_trisurf(dphif,detaf,pt3,cmap='jet', linewidths=0.5)

        plt.title(r'$p_{T}=3,Ridge$', size = 40)
        ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
        ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
        # axes.set_zlabel("$\Delta\eta$")
        plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
        plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

        plt.tight_layout()
        # plt.show()

        fig.savefig('pt=3_Ridge_3dplot.png')

        fig.clear()

    if numbers == 5:
        ax3d = plt.axes(projection="3d")
        pt4=np.loadtxt('pt_3dplot_Ridge.csv',delimiter=',',usecols=[6],skiprows=1)
        ax3d.plot_trisurf(dphif,detaf,pt4,cmap='jet', linewidths=0.5)

        plt.title(r'$p_{T}=4,Ridge$', size = 40)
        ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
        ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
        # axes.set_zlabel("$\Delta\eta$")
        plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
        plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

        plt.tight_layout()
        # plt.show()

        fig.savefig('pt=4_Ridge_3dplot.png')
        
if __name__ == '__main__':
    procs = []

    for i in range(1,11):
        proc = Process(target=functions, args=(i,))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()