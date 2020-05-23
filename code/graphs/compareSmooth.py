refFilePath = "../data/sapIngB1.fa"
kmerModelFilePath = "../data/kmer_model.hdf5"
sampleRead = "../data/pos-basecalled/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1706_ch_249_strand.fast5"

smoothParam = 5
repeatSignal = 10
basesSize = 50

import sys
import numpy as np
import matplotlib.pyplot as plt

import nadavca
from nadavca.dtw import KmerModel

sys.path.append("../helpers/hypothesis")
from signalHelper import (
    getSignalFromRead,
    smoothSignal,
    smoothSignalMed,
    stringToSignal,
)

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################
# nadavca
nadavca_align = nadavca.align_signal(
    refFilePath, [sampleRead], bwa_executable="./bwa/bwa"
)

assert len(nadavca_align) == 1, "Error! More than one alignment!"
nadavca_align = nadavca_align[0]

fromSignal, toSignal = nadavca_align[0].signal_range
table = nadavca_align[1][:basesSize]
refSeq = "".join(nadavca_align[0].reference_part)[:basesSize]

readSignal = getSignalFromRead(sampleRead)
readSignal = readSignal[table[0][2] : table[-1][2]]

refSignal = stringToSignal(refSeq, mod, repeatSignal=repeatSignal)

refSignal = np.array(refSignal, dtype=float)
readSignal = np.array(readSignal, dtype=float)

refSignal -= np.mean(refSignal)
readSignal -= np.mean(readSignal)

refSignal /= np.std(refSignal)
readSignal /= np.std(readSignal)

refSignalAvg = smoothSignal(refSignal, smoothParam)
readSignalAvg = smoothSignal(readSignal, smoothParam)

refSignalMed = smoothSignalMed(refSignal, smoothParam)
readSignalMed = smoothSignalMed(readSignal, smoothParam)

f, (ax1, ax2, ax3) = plt.subplots(3, 2)
ax1[0].plot(range(len(refSignal)), refSignal)
ax1[1].plot(range(len(readSignal)), readSignal)
ax2[0].plot(range(len(refSignalAvg)), refSignalAvg)
ax2[1].plot(range(len(readSignalAvg)), readSignalAvg)
ax3[0].plot(range(len(refSignalMed)), refSignalMed)
ax3[1].plot(range(len(readSignalMed)), readSignalMed)
plt.show()
