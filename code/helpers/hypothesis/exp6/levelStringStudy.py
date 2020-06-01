################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

readsPosFilePath = "../../../data/pos-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"
refFilePath = "../../../data/sapIngB1.fa"
alignedSquiggles = "../prepareData/alignedSquiggles2.txt"

smoothParam = 5
repeatSignal = 10
workingLen = 10000

readNum = 100

kmerLen = list(range(4, 36, 1))
levels = list(range(4, 15, 1))

plotLevels = [4, 5, 7, 9, 11, 13]
plotKmerLen = [4, 7, 13, 17, 21, 28]

import os
import sys
import math
import random
import h5py
import numpy as np
import mappy as mp
from pyfaidx import Fasta
import nadavca
from nadavca.dtw import KmerModel

import matplotlib
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

sys.path.append("../")
from signalHelper import (
    stringToSignal,
    getSignalFromRead,
    getReadsInFolder,
    stringAllignment,
    overlappingKmers,
    computeNorm,
    computeString,
    smoothSignal,
    countDashes,
    getAlignedIndex,
)


def plotAOC(src):
    src.sort()
    src.reverse()
    X, Y = [0], [0]
    x, y = 0, 0
    for i in src:
        if i[1] == 0:
            y += 1
        else:
            x += 1
        X.append(x)
        Y.append(y)
    return X, Y


def autolabel(ax, rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

################################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=0)

index = getAlignedIndex(alignedSquiggles)

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

pomerySm = [[] for i in levels]
pomery = [[] for i in levels]

readCounter = 0

for posRead in posReadsPaths:
    if readCounter == readNum:
        break

    readName = os.path.basename(posRead)
    
    if readName not in index:
        continue
    
    print(f"Working on {posRead}")
    
    fromSignal, toSignal = index[readName][4], index[readName][5]
    fromRef, toRef = index[readName][2], index[readName][3]
    strand = index[readName][1]
    ctg = index[readName][0]

    if (toSignal - fromSignal) < workingLen:
        continue

    #print(f"Signal alligned from {fromSignal} to {toSignal}")
    print("Working on", posRead)
    print(f"So far done {readCounter} reads")
    readCounter += 1

    if strand == 1:
        refSeq = str(Fasta(refFilePath)[ctg][fromRef:toRef])
    else:
        refSeq = str(-Fasta(refFilePath)[ctg][fromRef:toRef])

    refSignal = np.array(stringToSignal(refSeq, mod, repeatSignal=repeatSignal), float)
    readSignal = np.array(getSignalFromRead(posRead), dtype=float)
    readSignal = readSignal[fromSignal:toSignal]

    readSignal = readSignal[:workingLen]
    refSignal = refSignal[:workingLen]

    readSignalSm = smoothSignal(readSignal, smoothParam)
    refSignalSm = smoothSignal(refSignal, smoothParam)
    
    readShift, readScale = computeNorm(readSignal, 0, len(readSignal))
    readShiftSm, readScaleSm = computeNorm(readSignalSm, 0, len(readSignalSm))
    refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
    refShiftSm, refScaleSm = computeNorm(refSignalSm, 0, len(refSignalSm))
    
    readStrings, readStringsSm, refStrings, refStringsSm = {}, {}, {}, {}

    for l in levels:
        readStrings[l] = computeString(
            readSignal, 0, len(readSignal), readShift, readScale, l, overflow=0.30,
        )
        readStringsSm[l] = computeString(
            readSignalSm, 0, len(readSignalSm), readShiftSm, readScaleSm, l, overflow=0.30,
        )
        refStrings[l] = computeString(
            refSignal, 0, len(refSignal), refShift, refScale, l, overflow=0.30,
        )
        refStringsSm[l] = computeString(
            refSignalSm, 0, len(refSignalSm), refShiftSm, refScaleSm, l, overflow=0.30,
        )
        
        refLen, refSmLen = len(refStrings[l]), len(refStringsSm[l])
        readLen, readSmLen = len(readStrings[l]), len(readStringsSm[l])
        
        pomery[levels.index(l)].append(abs(refLen-readLen)/refLen)
        pomerySm[levels.index(l)].append(abs(refSmLen-readSmLen)/refSmLen)

dim1, dim2 = 2, len(plotLevels)//2
fig, axs = plt.subplots(dim1, dim2)

labels = [0, 5, 10, 30, 50]

for i in range(len(plotLevels)):
    X, Y = [], []
    data1 = pomerySm[levels.index(plotLevels[i])]
    data2 = pomery[levels.index(plotLevels[i])]
    
    good = [0] * len(labels)
    bad = [0] * len(labels)
    
    for x in data1:
        for k in reversed(range(len(labels))):
            if x*100 >= labels[k]:
                good[k] += 1
                break

    for x in data2:
        for k in reversed(range(len(labels))):
            if x*100 >= labels[k]:
                bad[k] += 1
                break
    
    #good = [sum([1 for x in data1 if x*100>labels[k]]) 
    #bad = [sum([1 for x in data2 if x*100>labels[k]]) for k in range(len(labels))]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    ax = axs[i//dim2, i%dim2] 

    rects1 = ax.bar(x - width/2, good, width, label='Smoothing')
    rects2 = ax.bar(x + width/2, bad, width, label='Not smoothing')

    if i%dim2 == 0:
        ax.set_ylabel('Frequency')
    if i//dim2 == 1:
        ax.set_xlabel("Percentage relative error (real to simulated)")
    ax.set_xticks(x)
    ax.set_xticklabels([">=" + str(k) + "%" for k in labels])
    #ax.legend()
    ax.set_ylim(top = readNum + 10)
    ax.set_title(f"{plotLevels[i]} levels")

    autolabel(ax, rects1)
    autolabel(ax, rects2)

#fig.suptitle(f'Length ratio of real squiggle level string to reference squiggle level string')

handles, labels = axs[dim1-1, dim2-1].get_legend_handles_labels()
fig.subplots_adjust(bottom=0.1, wspace=0.1)
fig.legend(handles, labels, loc='lower center', ncol=dim1*dim2)

plt.show()
