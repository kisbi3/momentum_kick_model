import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# mpl.rcParams["text.usetex"] = True

fig = plt.figure(figsize = (25,18))
ax = plt.axes()

fig.set_size_inches(23.384, 18, forward=True)

data = [0, 5.851025, 0.142829, 0.008921, 0.000728]
data_removejacob = [0, 5.883451, 0.14305, 0.008927, 0.000728]
data_remove_allE = [0, 0.59667, 0.071744, 0.005953, 0.000546]
data_remove_E = [0, 4.275051, 0.055528, 0.002325, 0.000143]
data_remove_Ei=[0, 0.816629, 0.194731, 0.024102, 0.00293]
data_remove_all=[0, 0.600328, 0.071855, 0.005957, 0.000546]
paper = [0, 1.23284673944207, 0.123284673944207, 0.010150641893705, 0.00094194593018]
dyratio = [1.411156, 0.694101, 0.413075, 0.325581, 0.277929]
#pt 0, 1, 2, 3, 4

print('pt\t\tratio\t\tratio/dyratio\tratio_remove_allE\tratio_remove_Ei')
for i in 1, 2, 3, 4:
    print(i,'\t',paper[i]/data[i], paper[i]/(data[i]*dyratio[i-1]), paper[i]/(data_remove_Ei[i]), paper[i]/(data_remove_Ei[i]*dyratio[i-1]))

print(data_remove_Ei[1]*dyratio[0], (data_remove_Ei[1]*dyratio[0]-paper[1])*100/paper[1])
print(data_remove_Ei[2]*dyratio[1], (data_remove_Ei[2]*dyratio[1]-paper[2])*100/paper[2])
print(data_remove_Ei[3]*dyratio[2], (data_remove_Ei[3]*dyratio[2]-paper[3])*100/paper[3])
print(data_remove_Ei[4]*dyratio[3], (data_remove_Ei[4]*dyratio[3]-paper[4])*100/paper[4])


