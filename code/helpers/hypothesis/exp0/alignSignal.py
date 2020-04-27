################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

import sys

# here we want to chose some part of signal that is not close to beg or end of read
refFilePath = "../../data/sapIngB1.fa"
sampleRead = sys.argv[1]
#sampleFakeRead = args.sampleFakeRead

# here we want to chose some part of signal that is not close to beg or end of read
fromSignal, toSignal = 20000, 22000
#fromSignalFake, toSignalFake = 15000, 17000

kmerModelFilePath = "../../data/kmer_model.hdf5"
# number of levels we use
numLevels = 6
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 8
kernelLen = 10
winSize = 15 * kernelLen

import math
import h5py
import copy
import numpy as np
import mappy as mp
import nadavca
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import getLevels, getLevelString, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead, Table_Iterator, seqSignalCor

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

############s####################################################################
# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)
nadavca_align = nadavca.align_signal(refFilePath, [sampleRead], bwa_executable='./bwa/bwa')

assert (len(nadavca_align) == 1), "Error! More than one alignment!"

nadavca_align = nadavca_align[0]

assert (nadavca_align[0].reverse_complement == False), "Error! Reverse strand!"

fromSignal, toSignal = nadavca_align[0].signal_range

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)
originalSignal = originalSignal[fromSignal:toSignal]

# create artificial signal
signalFrTo = "".join(nadavca_align[0].reference_part)
artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)

# load some fake signal
#fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = copy.deepcopy(originalSignal[::-1])
fakeSignal = np.array(fakeSignal, dtype = float)
################################################################################
# get levels

normalizeWindow(originalSignal)
normalizeWindow(artifSignal)
normalizeWindow(fakeSignal)

levelA = getLevels(artifSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)

# count have many strings are shared by fabricatedSignal and original/processedSignal
levelO = getLevels(originalSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
counterO = len([i for i in levelO if i in levelA])
print("Artificial match rate with original signal {0}/{1}".format(counterO, len(levelO)))

levelF = getLevels(fakeSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
counterF = len([i for i in levelF if i in levelA])
print("Total fake match rate with original signal {0}/{1}".format(counterF, len(levelF)))

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

################################################################################
# ploting

from matplotlib.widgets import Slider, Button, TextBox

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
    origW = copy.deepcopy(originalSignal[newBeg1:(newBeg1+winSize)])
    helperArr = nadavca_align[0].alignment
    
    best, bestB = 1000000, -1
    for i in nadavca_align[0].alignment:
        if i[0] > sOrig.val:
            best = i[1]
            break
    
    print(best*repeatSignal)
    artifW = copy.deepcopy(artifSignal[newBeg2:(newBeg2+winSize)])
    shiftOrigUp = copy.deepcopy(origW)
    shiftOrigDw = copy.deepcopy(origW)
    normalizeWindow(origW)
    normalizeWindow(artifW)
    normalizeWindow(shiftOrigUp, shift = 0.4)
    normalizeWindow(shiftOrigDw, shift = -0.4)
    l1.set_ydata(origW)
    l2.set_ydata(artifW)
    levelS1 = getLevelString(origW, numLevels = numLevels)
    levelS2 = getLevelString(artifW, numLevels = numLevels)
    levelS3 = getLevelString(shiftOrigUp, numLevels = numLevels)
    levelS4 = getLevelString(shiftOrigDw, numLevels = numLevels)
    # orig
    print("O:" + levelS1)
    # artif
    print("A:" + levelS2)
    # shift up
    print("U:" + levelS3)
    # shift down
    print("D:" + levelS4)
    fig.canvas.draw_idle()

sOrig.on_changed(update)
sArtif.on_changed(update)

def daco(val):
    print("Ehe")
    sArtif.val = float(val)
    update(val)

axbox = plt.axes([0.7, 0.05, 0.1, 0.075])
text_box = TextBox(axbox, 'Evaluate', initial="0.0")
text_box.on_submit(daco)

plt.show()
