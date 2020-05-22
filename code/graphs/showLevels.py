refFilePath = "../data/sapIngB1.fa"
read = "../data/pos-basecalled/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1706_ch_249_strand.fast5"

#fromSignal, toSignal = 10050, 10110
fromSignal, toSignal = 10050, 10200
#fromSignal, toSignal = 10150, 10200
levels = 4

import sys
import numpy as np

sys.path.append("../helpers/hypothesis")
from signalHelper import getSignalFromRead

import matplotlib.pyplot as plt

signal = getSignalFromRead(read)[fromSignal:toSignal]

mini, maxi = min(signal), max(signal)
levelSize = (maxi-mini)/levels

y_values = [chr(ord("a") + i) for i in range(4)]
y_axis = np.arange(0, 4, 1)

for a in np.arange(mini, maxi+0.00001, levelSize):
    plt.axhline(y=a, color = 'r', linewidth = '2')

helper = [chr(ord("a") + int((i-mini)/levelSize)) for i in signal]
signalLevels = [" "] * len(helper)

events = []
begg = 0
for i in range(len(helper)):
    if helper[i] != helper[begg]:
        events.append((begg,i))
        begg = i
events.append((begg, i))

for i in range(len(events)):
    event = events[i]
    signalLevels[(event[0]+event[1])//2] = helper[event[0]]
    c = ["#e28c06", "#e10304", "#b61aae", "#dbfd45"][ord(helper[event[0]])-ord("a")-1]
    #plt.fill_between(range(event[0], event[1]), signal[event[0]:event[1]], y2 = mini, color = c)

plt.ylim(mini-0.1, maxi+0.1)

plt.yticks(np.arange(mini+levelSize/2, maxi+levelSize/2, levelSize), y_values)
plt.xticks(range(len(signal)), signalLevels)

for i in range(len(events)):
    event = events[i]
    c=np.random.rand(3,)
    b = event[0]
    e = event[1]
    plt.plot(range(e-1, e+1), signal[e-1:e+1], color = "#000000", linewidth="1", zorder=2)
    c = ["#e28c06", "#00ff00", "#00bfff", "#dbfd45"]
    plt.plot(range(b, e), signal[b:e], color = c[ord(helper[b])-ord("a")-1], linewidth = "3.0", zorder=1)
    #print(x)
    #plt.axvline(x=event[0], ymin = 0, ymax = b)
    #plt.axvline(x=event[1]-1, ymin = 0, ymax = e)

plt.scatter(range(len(signal)), signal, marker = 'o', color = "#000000", zorder=3)
plt.show()
