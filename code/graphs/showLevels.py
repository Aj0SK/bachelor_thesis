refFilePath = "../data/sapIngB1.fa"
read = "../data/pos-basecalled/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1706_ch_249_strand.fast5"

# fromSignal, toSignal = 10050, 10110
# fromSignal, toSignal = 10050, 10200
fromSignal, toSignal = 10150, 10200
levels = 6

import sys
import numpy as np

sys.path.append("../helpers/hypothesis")
from signalHelper import getSignalFromRead

import matplotlib.pyplot as plt

signal = getSignalFromRead(read)[fromSignal:toSignal]

mini, maxi = min(signal), max(signal) + 5
levelSize = (maxi - mini) / levels

y_values = [chr(ord("a") + i) for i in range(levels)]
y_axis = np.arange(0, levels, 1)

for a in np.arange(mini, maxi + levelSize, levelSize):
    plt.axhline(y=a, color="r", linewidth="2")

helper = [chr(ord("a") + int((i - mini) / levelSize)) for i in signal]
signalLevels = [" "] * len(helper)

events = []
begg = 0
for i in range(len(helper)):
    if helper[i] != helper[begg]:
        events.append((begg, i))
        begg = i
events.append((begg, i))

for i in range(len(events)):
    event = events[i]
    signalLevels[(event[0] + event[1]) // 2] = helper[event[0]]

plt.ylim(mini - 5, maxi + 5)

plt.yticks(np.arange(mini + levelSize / 2, maxi + levelSize / 2, levelSize), y_values)
plt.xticks(range(len(signal)), signalLevels)

c = ["#e28c06", "#00ff00", "#00bfff", "#dbfd45", np.random.rand(3,), np.random.rand(3,)]

for i in range(len(events)):
    event = events[i]
    b = event[0]
    e = event[1]
    plt.plot(
        range(e - 1, e + 1),
        signal[e - 1 : e + 1],
        color="#000000",
        linewidth="1",
        zorder=2,
    )
    plt.plot(range(b, e), signal[b:e], color="#000000", linewidth="1.0", zorder=1)
    plt.scatter(
        range(b, e),
        signal[b:e],
        marker="o",
        color=c[ord(helper[b]) - ord("a") - 1],
        zorder=3,
    )
    # print(x)
    # plt.axvline(x=event[0], ymin = 0, ymax = b)
    # plt.axvline(x=event[1]-1, ymin = 0, ymax = e)

plt.xlabel("Level string", fontsize=16)
plt.ylabel("Levels", fontsize=16)
# plt.scatter(range(len(signal)), signal, marker = 'o', color = "#000000", zorder=3)
plt.show()
