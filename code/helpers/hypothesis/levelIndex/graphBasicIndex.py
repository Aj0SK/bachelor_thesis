dataTable = {}
levels = []
kmerLens = []

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

def plotAOC(src):
    src.sort()
    src.reverse()
    X, Y = [0], [0]
    x, y = 0, 0
    for i in src:
        if i[1] == 1:
            y += 1
        else:
            x += 1
        X.append(x)
        Y.append(y)
    return X, Y

with open("imp.txt", "r") as f:
    for line in f.readlines():
        line = [i.strip() for i in line.split()]
        
        
        if line[0] != "Levels:":
            print("Warn")
            continue

        level = int(line[1])
        kmerLen = int(line[3])
        
        if line[4] == "#":
            data = [int(i) for i in line[5:]]
            data = [(data[i], data[i+1]) for i in range(0, len(data), 2)]
        
        dataTable[(level, kmerLen)] = data
        print(f"{level} {kmerLen}")
        
        levels.append(level)
        kmerLens.append(kmerLen)

levels = list(set(levels))
kmerLens = list(set(kmerLens))

dim1, dim2 = 1, 3
fig, axs = plt.subplots(dim1, dim2)
counter = 0

for l in levels:
    ax = axs[counter]
    for k in kmerLens:
        if (l, k) not in dataTable:
            continue
        data = dataTable[(l, k)]
        print(data)
        print(sum(1 for i in data if i[1] == 1))
        x, y = plotAOC(data)
        ax.plot(x, y, linewidth = 2, label=str(k))
        ax.plot(x, x, linestyle = '--')
    counter += 1

handles, labels = axs[dim2-1].get_legend_handles_labels()
fig.subplots_adjust(bottom=0.1, wspace=0.1)
leg = fig.legend(handles, labels, loc='lower center', ncol=dim1*dim2)

for line in leg.get_lines():
    line.set_linewidth(4.0)

plt.show()
