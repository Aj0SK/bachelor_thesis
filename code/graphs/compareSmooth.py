refFilePath = "../data/sapIngB1.fa"
read = "../data/pos-basecalled/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1706_ch_249_strand.fast5"

fromSignal, toSignal = 10050, 10200
smoothParam = 5

import sys
import numpy as np

sys.path.append("../helpers/hypothesis")
from signalHelper import getSignalFromRead, smoothSignal

import matplotlib.pyplot as plt

signal = getSignalFromRead(read)[fromSignal:toSignal]

signalAvg1 = smoothSignal(signal, smoothParam)

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(signal)
ax2.plot(signalAvg1)
plt.show()
