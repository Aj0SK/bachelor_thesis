################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

readsPosFilePath = "../../data/pos-basecalled"
readsNegFilePath = "../../data/neg-basecalled"
kmerModelFilePath = "../../data/kmer_model.hdf5"
refFilePath = "../../data/sapIngB1.fa"

repeatSignal = 10

levels = 12
minSignal, maxSignal = -2.0, 2.0

import math
import h5py
import numpy as np
import mappy as mp
from pyfaidx import Fasta
import nadavca
from nadavca.dtw import KmerModel

from signalHelper import stringToSignal, getSignalFromRead, getReadsInFolder, stringAllignment
from signalHelper import computeNorm, computeString, smoothSignal

import matplotlib
import matplotlib.pyplot as plt

############s###################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=0)
negReadsPaths = getReadsInFolder(readsNegFilePath, minSize=0)

sampleRead = posReadsPaths[0]
sampleFakeRead = negReadsPaths[0]

############s###################################################################
# load reference sequence and create index for fast mapping
ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"
################################################################################
# nadavca
nadavca_align = nadavca.align_signal(refFilePath, [sampleRead],
                                     bwa_executable='./bwa/bwa')

assert (len(nadavca_align) == 1), "Error! More than one alignment!"
nadavca_align = nadavca_align[0]
assert (nadavca_align[0].reverse_complement == False), "Error! Reverse strand!"
################################################################################

refStr = str(ref[nadavca_align[0].contig_name])
fromSignal, toSignal = nadavca_align[0].signal_range

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype=float)

table = nadavca_align[1][:40]
refSeq = "".join(nadavca_align[0].reference_part)[:40]

x, y = [], []

for entry in table:  #entry is list of [ref_index, signal_start, signal_end]
    x.append(str(refStr[entry[0]]))
    for i in range(entry[1], entry[2]):
        x.append(" ")
        y.append(originalSignal[i])

#plt.plot(y)
#plt.xticks(y_pos, x, color='orange', rotation=45, fontweight='bold', horizontalalignment='right')
#plt.tick_params(labelbottom='off')

refSignal = stringToSignal(refSeq, mod, repeatSignal=repeatSignal)
refSeqHelper = []

for i in refSeq:
    refSeqHelper.append(i)
    for k in range(repeatSignal - 1):
        refSeqHelper.append("_")

y = smoothSignal(y, 5)
refSignal = smoothSignal(refSignal, 5)

ySignalShift, ySignalScale = computeNorm(y, 0, len(y))
#y -= ySignalShift
#y /= ySignalScale

refSignalShift, refSignalScale = computeNorm(refSignal, 0, len(refSignal))
#refSignal -= refSignalShift
#refSignal /= refSignalScale

readString = computeString(y,
                           0,
                           len(y),
                           ySignalShift,
                           ySignalScale,
                           levels,
                           overflow=0.30)
refString = computeString(refSignal,
                          0,
                          len(refSignal),
                          refSignalShift,
                          refSignalScale,
                          levels,
                          overflow=0.30)

fakeSignal = getSignalFromRead(sampleFakeRead)[6000:6000+len(y)]
fakeSignal = smoothSignal(fakeSignal, 5)
fakeSignalShift, fakeSignalScale = computeNorm(fakeSignal, 0, len(fakeSignal))

fakeString = computeString(fakeSignal,
                          0,
                          len(fakeSignal),
                          fakeSignalShift,
                          fakeSignalScale,
                          levels,
                          overflow=0.30)

print("refstring vs Readstring vs fakestring")
print(refString)
print(readString)
print(fakeString)

print("Readstring allignment")
a, b = stringAllignment(refString, readString)
print(a)
print(b)

print("Fakestring allignment")
c, d = stringAllignment(refString, fakeString)
print(c)
print(d)

y = np.array(y, float)
refSignal = np.array(refSignal, float)

y -= ySignalShift
y /= ySignalScale

y[y<minSignal] = minSignal
y[y>maxSignal] = maxSignal

refSignal[refSignal<minSignal] = minSignal
refSignal[refSignal>maxSignal] = maxSignal

fig, axs = plt.subplots(2)
fig.suptitle('Read vs reference')
axs[0].plot(y)
axs[0].set_xticks(np.arange(len(x)))
axs[0].set_xticklabels(x)

axs[1].plot(refSignal)
axs[1].set_xticks(np.arange(len(refSeqHelper)))
axs[1].set_xticklabels(refSeqHelper)

for i in np.arange(minSignal, maxSignal + 0.001,
                   (maxSignal - minSignal) / levels):
    axs[0].axhline(y=i)
    axs[1].axhline(y=i)

plt.show()
