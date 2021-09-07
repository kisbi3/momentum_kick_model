import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import csv





f = open('pardis.csv', 'r', encoding='utf-8')
lines=[]
rdr = csv.reader(f, delimiter = "\t")
# print(lines)
for line in rdr:
    lines.append(line)
# print(lines)

x=[]
y=[]

for line in lines[4:93]:
    # x.append(float(line[0]))
    # y.append(float(line[1]))
    x.append(line[0].split(' '))
    # y.append(line[1])
    # Ncoll.append(float(line[2]))
print(x)


f.close()