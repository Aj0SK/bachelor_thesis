################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

# here we want to chose some part of signal that is not close to beg or end of read
fromSignal, toSignal = 19000, 21000
fromSignalFake, toSignalFake = 15000, 17000
refFilePath = "../../data/sapIngB1.fa"
sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"
sampleFakeRead = "../../data/albacore-output-neg/magnu_20181218_FAH93149_MN26672_sequencing_run_sapFun_DNase_flush_88184_read_1001_ch_42_strand.fast5"

kmerModelFilePath = "../../data/kmer_model.hdf5"
# number of levels we use
numLevels = 6
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 9
kernelLen = 7
winSize = 15 * kernelLen

import math
import h5py
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import getLevels, getLevelString, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead, Table_Iterator

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

################################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

################################################################################
# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastq, basecallTable = getSeqfromRead(sampleRead)
# read second line and decode it into string from bytes
fastqSeq = fastq.splitlines()[1].decode()

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

################################################################################
# get only part of signal and find corresponding basecalled sequence, then find
# it in reference
signalFrTo = ""

for i in Table_Iterator(basecallTable):
    #print(chr(i[0]) + " " + str(i[1]) + " " + str(i[2]))
    if i[1] >= fromSignal and i[1] <= toSignal:
        signalFrTo += str(chr(i[0]))

#print(signalFrTo)

# try to find string in reference
hits = list(sequenceIndex.map(signalFrTo))
assert len(hits) == 1, "Too many hits"
signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])

# create artificial signal
artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)
originalSignal = originalSignal[fromSignal:toSignal]

################################################################################
# modify original signal
from scipy.signal import medfilt

procSignal = np.array(getSignalFromRead(sampleRead), dtype = float)[fromSignal:toSignal]
procSignal = medfilt(procSignal, kernel_size = 9)

levelO = getLevels(originalSignal, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
levelO = set(levelO)
levelE = getLevels(procSignal, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
levelE = set(levelE)

levelA = getLevels(artifSignal, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)

# count have many strings are shared by fabricatedSignal and original/processedSignal
counterO = len([i for i in levelA if i in levelO])
counterE = len([i for i in levelA if i in levelE])

print("Artificial match rate with original signal {0}/{1}".format(counterO, len(levelA)))
print("Artificial match rate with preprocessed signal {0}/{1}".format(counterE, len(levelA)))

# load original signal from read
fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = np.array(fakeSignal, dtype = float)
fakeSignal = fakeSignal[fromSignalFake:toSignalFake]
levelF = getLevels(fakeSignal, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
counterF = len([i for i in levelA if i in levelF])
print("Total fake match rate with original signal {0}/{1}".format(counterF, len(levelA)))

#print(levelO)
#print(levelE)
#print(levelF)

################################################################################
#

counter = 0
for w in levelO:
    for i in range(kernelLen-1):
        if abs(ord(w[i])-ord(w[i+1])) > 1:
            counter += 1

print("Jumps {0}/{1}".format(counter, len(levelO)*kernelLen))

counter = 0
for w in levelF:
    for i in range(kernelLen-1):
        if abs(ord(w[i])-ord(w[i+1])) > 1:
            counter += 1

print("Jumps {0}/{1}".format(counter, len(levelF)*kernelLen))
################################################################################
# Graph normalized original, processed and artificial signal

#normalizeWindow(originalSignal)
#normalizeWindow(artifSignal)
#normalizeWindow(procSignal)
#normalizeWindow(fakeSignal)

o = np.linspace(0, originalSignal.shape[0], originalSignal.shape[0])
a = np.linspace(0, artifSignal.shape[0], artifSignal.shape[0])
p = np.linspace(0, procSignal.shape[0], procSignal.shape[0])
f = np.linspace(0, fakeSignal.shape[0], fakeSignal.shape[0])

fig, axs = plt.subplots(3)
fig.suptitle('Original, modified original(processed) and artificial signal')
axs[0].plot(o, originalSignal, label='data 1', color = 'r')
#axs[1].plot(p, procSignal, label='data 2', color = 'g')
axs[1].plot(f, fakeSignal, label='data 2', color = 'g')
axs[2].plot(a, artifSignal, label='data 3', color = 'b')

axs[0].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(o), linewidth=1, color='gray')
axs[1].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(p), linewidth=1, color='gray')
axs[2].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(a), linewidth=1, color='gray')
plt.show(block=False)

from matplotlib.widgets import Slider

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
plt.subplots_adjust(left=0.25, bottom=0.25)
l1, = ax1.plot(np.arange(winSize), originalSignal[:winSize], lw=2)
l2, = ax2.plot(np.arange(winSize), artifSignal[:winSize], lw=2)
#ax1.margins(x=0)

ax1.set_ylim([minSignal-0.1, maxSignal+0.1])
ax2.set_ylim([minSignal-0.1, maxSignal+0.1])

axcolor = 'lightgoldenrodyellow'
axorig = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axartif = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sOrig = Slider(axorig, 'Orig', 0, toSignal-fromSignal-winSize, valinit=0, valstep=1)
sArtif = Slider(axartif, 'Artif', 0, artifSignal.shape[0]-winSize, valinit=0, valstep=1)

def update(val):
    newBeg1, newBeg2 = int(sOrig.val), int(sArtif.val)
    w1 = copy.deepcopy(originalSignal[newBeg1:(newBeg1+winSize)])
    w2 = copy.deepcopy(artifSignal[newBeg2:(newBeg2+winSize)])
    normalizeWindow(w1)
    normalizeWindow(w2)
    #l1.set_xdata(np.arange(newBeg1, newBeg1+winSize))
    #l2.set_xdata(np.arange(newBeg2, newBeg2+winSize))
    l1.set_ydata(w1)
    l2.set_ydata(w2)
    print(getLevelString(w1))
    print(getLevelString(w2))
    fig.canvas.draw_idle()

sOrig.on_changed(update)
sArtif.on_changed(update)

plt.show()
