################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

# here we want to chose some part of signal that is not close to beg or end of read
fromSignal, toSignal = 20000, 24000
fromSignalFake, toSignalFake = 15000, 19000
refFilePath = "../../data/sapIngB1.fa"
sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1007_ch_424_strand.fast5"
#sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"
sampleFakeRead = "../../data/albacore-output-neg/magnu_20181218_FAH93149_MN26672_sequencing_run_sapFun_DNase_flush_88184_read_1001_ch_42_strand.fast5"

kmerModelFilePath = "../../data/kmer_model.hdf5"
# number of levels we use
numLevels = 4
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 8
kernelLen = 7
winSize = 15 * kernelLen

import math
import h5py
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import getLevels, getLevelString, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead, Table_Iterator, seqSignalCor

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

############s####################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

################################################################################
# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastqSeq, basecallTable = getSeqfromRead(sampleRead)
# read second line and decode it into string from bytes

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

################################################################################
# get only part of signal and find corresponding basecalled sequence, then find
# it in reference
signalFrTo = seqSignalCor(fromSignal, toSignal, basecallTable)

#print(signalFrTo)

# try to find string in reference
hits = list(sequenceIndex.map(signalFrTo))
assert len(hits) == 1, "Too many hits"
signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])

# create artificial signal
artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)
originalSignal = originalSignal[fromSignal:toSignal]

# load some fake signal
fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = np.array(fakeSignal, dtype = float)
fakeSignal = fakeSignal[fromSignalFake:toSignalFake]

################################################################################
# modify original signal
from scipy.signal import medfilt

levelO = getLevels(originalSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
levelO = set(levelO)

levelA = getLevels(artifSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)

# count have many strings are shared by fabricatedSignal and original/processedSignal
counterO = len([i for i in levelA if i in levelO])

print("Artificial match rate with original signal {0}/{1}".format(counterO, len(levelA)))

# load original signal from read
fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = np.array(fakeSignal, dtype = float)
fakeSignal = fakeSignal[fromSignalFake:toSignalFake]
levelF = getLevels(fakeSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
counterF = len([i for i in levelA if i in levelF])
print("Total fake match rate with original signal {0}/{1}".format(counterF, len(levelA)))

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
#

'''
levelStringStarts = {}
for winBeg in range(0, originalSignal.shape[0]-winSize, 100):
        winEnd = winBeg + winSize
        currWindow = copy.deepcopy(originalSignal[winBeg:winEnd])
        normalizeWindow(currWindow, minSignal, maxSignal)
        outString = getLevelString(currWindow, minSignal, maxSignal, numLevels)
        if len(outString) < kernelLen:
            continue
        outString = outString[:kernelLen]
        levelStringStarts[outString] = winBeg

fastq1, basecallTable1 = getSeqfromRead(sampleRead)
fastq2, basecallTable2 = getSeqfromRead(sampleFakeRead)

for winBeg in range(0, artifSignal.shape[0]-winSize, 1):
        winEnd = winBeg + winSize

        # normalize window
        currWindow = copy.deepcopy(artifSignal[winBeg:winEnd])
        normalizeWindow(currWindow, minSignal, maxSignal)

        # cut into horizontal levels
        outString = getLevelString(currWindow, minSignal, maxSignal, numLevels)
        if len(outString) < kernelLen:
            continue
        outString = outString[:kernelLen]
        if outString in levelStringStarts:
            print("Zhoda: {0} {1}".format(winBeg, levelStringStarts[outString]))
            print(outString)
            #f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
            sig1 = copy.deepcopy(originalSignal[levelStringStarts[outString]:levelStringStarts[outString]+winSize])
            sig2 = copy.deepcopy(artifSignal[winBeg:winEnd])
            normalizeWindow(sig1)
            normalizeWindow(sig2)
            #l1, = ax1.plot(np.arange(winSize), sig1, lw=2)
            #l2, = ax2.plot(np.arange(winSize), sig2, lw=2)
            #ax1.set_ylim([minSignal-0.1, maxSignal+0.1])
            #ax2.set_ylim([minSignal-0.1, maxSignal+0.1])
            #ax1.hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=winSize, linewidth=1, color='gray')
            #ax2.hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=winSize, linewidth=1, color='gray')
            #plt.show()
            signalFrTo1 = seqSignalCor(levelStringStarts[outString], levelStringStarts[outString]+winSize, basecallTable1)
            signalFrTo2 = seqSignalCor(winBeg, winEnd, basecallTable2)
            print(signalFrTo1)
            print(signalFrTo2)

'''
        

################################################################################
# Graph normalized original, processed and artificial signal

#normalizeWindow(originalSignal)
#normalizeWindow(artifSignal)
#normalizeWindow(fakeSignal)

o = np.linspace(0, originalSignal.shape[0], originalSignal.shape[0])
a = np.linspace(0, artifSignal.shape[0], artifSignal.shape[0])
f = np.linspace(0, fakeSignal.shape[0], fakeSignal.shape[0])

fig, axs = plt.subplots(3)
fig.suptitle('Original, modified original(processed) and artificial signal')
axs[0].plot(o, originalSignal, label='data 1', color = 'r')
axs[1].plot(f, fakeSignal, label='data 2', color = 'g')
axs[2].plot(a, artifSignal, label='data 3', color = 'b')

axs[0].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(o), linewidth=1, color='gray')
axs[2].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(a), linewidth=1, color='gray')
plt.show(block=False)

from matplotlib.widgets import Slider
import edlib

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
plt.subplots_adjust(left=0.25, bottom=0.25)
l1, = ax1.plot(np.arange(winSize), originalSignal[:winSize], lw=2)
l2, = ax2.plot(np.arange(winSize), artifSignal[:winSize], lw=2)

ax1.set_ylim([minSignal-0.1, maxSignal+0.1])
ax2.set_ylim([minSignal-0.1, maxSignal+0.1])

ax1.hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=winSize, linewidth=1, color='gray')
ax2.hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=winSize, linewidth=1, color='gray')

axcolor = 'lightgoldenrodyellow'
axorig = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axartif = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sOrig = Slider(axorig, 'Orig', 0, toSignal-fromSignal-winSize, valinit=0, valstep=1)
sArtif = Slider(axartif, 'Artif', 0, artifSignal.shape[0]-winSize, valinit=0, valstep=1)

def update(val):
    newBeg1, newBeg2 = int(sOrig.val), int(sArtif.val)
    w1 = copy.deepcopy(originalSignal[newBeg1:(newBeg1+winSize)])
    w2 = copy.deepcopy(artifSignal[newBeg2:(newBeg2+winSize)])
    w3 = copy.deepcopy(w1)
    w4 = copy.deepcopy(w1)
    normalizeWindow(w1)
    normalizeWindow(w2)
    normalizeWindow(w3, shift = 0.4)
    normalizeWindow(w4, shift = -0.4)
    l1.set_ydata(w1)
    l2.set_ydata(w2)
    levelS1 = getLevelString(w1, numLevels = numLevels)
    levelS2 = getLevelString(w2, numLevels = numLevels)
    levelS3 = getLevelString(w3, numLevels = numLevels)
    levelS4 = getLevelString(w4, numLevels = numLevels)
    #if len(levelS1) >= kernelLen:
    #    levelS1 = levelS1[:kernelLen]
    #if len(levelS2) >= kernelLen:
    #    levelS2 = levelS2[:kernelLen]
    #
    print(levelS1)
    print(levelS2)
    print(levelS3)
    print(levelS4)
    print("1. distance is {0}".format(edlib.align(levelS1, levelS2)["editDistance"]))
    print("2. distance is {0}".format(edlib.align(levelS3, levelS2)["editDistance"]))
    print("3. distance is {0}".format(edlib.align(levelS4, levelS2)["editDistance"]))
    fig.canvas.draw_idle()

sOrig.on_changed(update)
sArtif.on_changed(update)

plt.show()
