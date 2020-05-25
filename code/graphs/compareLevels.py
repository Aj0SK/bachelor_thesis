refFilePath = "../data/sapIngB1.fa"
read = "../data/pos-basecalled/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1706_ch_249_strand.fast5"

workLen = 50
levels = 3
mini, maxi = -2.0, 2.0

import sys
import numpy as np
from itertools import groupby

sys.path.append("../helpers/hypothesis")
from signalHelper import getSignalFromRead, stringAllignment

import matplotlib.pyplot as plt

def normal(signal):
    newsignal = signal - np.mean(signal)
    newsignal /= np.std(newsignal)
    return newsignal

def f(signal, mini, levelSize):
    helper = [chr(ord("a") + int((i-mini)/levelSize)) for i in signal]
    helper = "".join(helper)
    helper = "".join([k for k, g in groupby(helper)])
    return helper

signal = getSignalFromRead(read)
signal1, signal2 = [], []
found = False

levelSize = 0, 0, 0

signal = normal(signal)

signal[signal>maxi] = maxi
signal[signal<mini] = mini

helper = {}

counter = 3

for i in range(0, len(signal)-2*workLen+1, 30):
    signal1 = signal[i:i+workLen]

    levelSize = (maxi-mini)/levels

    str1 = f(signal1, mini, levelSize)

    if str1 in helper:
        counter -= 1
        if counter != 0:
            continue
        print(i)
        print(helper[str1])
        signal2 = signal[helper[str1]:helper[str1]+workLen]
        #signal2 = normal(signal2)
        break
    
    helper[str1] = i

#signal2[20:] += 0.7
#signal1[:20] -= 0.5

for a in np.arange(mini, maxi+levelSize, levelSize):
    plt.axhline(y=a, color = 'r', linewidth = '2')
plt.plot(signal1)
plt.plot(signal2)
plt.show()

import nadavca
from nadavca.dtw import KmerModel

nadavca_align = nadavca.align_signal(
    refFilePath, [read], bwa_executable="./bwa/bwa"
)

assert len(nadavca_align) == 1, "Error! More than one alignment!"
nadavca_align = nadavca_align[0]

fromSignal, toSignal = nadavca_align[0].signal_range
table = nadavca_align[1][:30]
refSeq = "".join(nadavca_align[0].reference_part)[:30]

signal = signal[table[0][1]:table[10][2]]
signal = normal(signal)

signal[signal<mini] = mini
signal[signal>maxi] = maxi

for a in np.arange(mini, maxi, (maxi-mini)/12):
    plt.axhline(y=a, color = 'r', linewidth = '2')
plt.plot(signal)
plt.show()
