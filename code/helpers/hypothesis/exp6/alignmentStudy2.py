################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"
refFilePath = "../../../data/sapIngB1.fa"
alignedSquiggles = "../prepareData/alignedSquiggles.txt"

smoothParam = 5
repeatSignal = 10
workingLen = 5000

readNum = 200

kmerLen = list(range(4, 36, 1))
levels = list(range(4, 15, 1))

plotLevels = [4, 5, 7, 9, 11, 13]  # range(4, 10)
plotKmerLen = [4, 7, 13, 17, 21, 28]  # range(4, 40, 10)

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
    getAlignedIndex
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


################################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=0)
negReadsPaths = getReadsInFolder(readsNegFilePath, minSize=0)

index = getAlignedIndex(alignedSquiggles)

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

pomery = [[[] for j in kmerLen] for i in levels]
overlap = [[[] for j in kmerLen] for i in levels]
goodDash = [[] for i in levels]
badDash = [[] for i in levels]
alignLenRead = [0 for _ in levels]
alignLenFake = [0 for _ in levels]
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
    fakeSignal = []
    fakeIndex = -1
    while len(fakeSignal) <= toSignal:
        fakeIndex = random.randint(0, len(negReadsPaths) - 1)
        fakeSignal = np.array(getSignalFromRead(negReadsPaths[fakeIndex]), dtype=float)
    fakeSignal = fakeSignal[fromSignal:toSignal]

    readSignal = readSignal[:workingLen]
    refSignal = refSignal[:workingLen]
    fakeSignal = fakeSignal[:workingLen]

    readSignal = smoothSignal(readSignal, smoothParam)
    refSignal = smoothSignal(refSignal, smoothParam)
    fakeSignal = smoothSignal(fakeSignal, smoothParam)
    readShift, readScale = computeNorm(readSignal, 0, len(readSignal))
    refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
    fakeShift, fakeScale = computeNorm(fakeSignal, 0, len(fakeSignal))

    readStrings, refStrings, fakeStrings = {}, {}, {}

    for l in levels:
        readStrings[l] = computeString(
            readSignal, 0, len(readSignal), readShift, readScale, l, overflow=0.30,
        )
        refStrings[l] = computeString(
            refSignal, 0, len(refSignal), refShift, refScale, l, overflow=0.30,
        )
        fakeStrings[l] = computeString(
            fakeSignal, 0, len(fakeSignal), fakeShift, fakeScale, l, overflow=0.30,
        )

    for l in levels:
        a, b = stringAllignment(refStrings[l], readStrings[l])
        c, d = stringAllignment(refStrings[l], fakeStrings[l])
        
        #alignLenRead[levels.index(l)] += len(a)
        #alignLenFake[levels.index(l)] += len(c) 
        
        alignLenRead[levels.index(l)] += len(readStrings[l])
        alignLenFake[levels.index(l)] += len(fakeStrings[l]) 
        
        #a = a[:300]
        #b = b[:300]
        #c = c[:300]
        #d = d[:300]
        
        dashes1 = [countDashes(a, i) + countDashes(b, i) for i in range(1, 21)]
        dashes2 = [countDashes(c, i) + countDashes(d, i) for i in range(1, 21)]
        
        goodDash[levels.index(l)].append(dashes1)
        badDash[levels.index(l)].append(dashes2)
        
        #print(dashes1)
        #print(dashes2)
        
        '''
        for k in kmerLen:
            counter = 0
            for i in range(len(readStrings[l]) - k + 1):
                for j in range(len(refStrings[l]) - k + 1):
                    w1 = readStrings[l][i : i + k]
                    w2 = refStrings[l][j : j + k]
                    if ("-" not in w1) and ("-" not in w2):
                        counter += 1'''

    for l in levels:
        print(f"levels({l}): #", end="")
        for k in kmerLen:
            goodOverlap = overlappingKmers(refStrings[l], readStrings[l], k)
            badOverlap = overlappingKmers(refStrings[l], fakeStrings[l], k)

            overlap[levels.index(l)][kmerLen.index(k)].append((goodOverlap, 0))
            overlap[levels.index(l)][kmerLen.index(k)].append((badOverlap, 1))

            if goodOverlap == 0:
                badOverlap = 1000
                goodOverlap = 100
            if (badOverlap / goodOverlap) > 2.0:
                badOverlap = 2 * goodOverlap
            pomery[levels.index(l)][kmerLen.index(k)].append(
                100.0 * badOverlap / goodOverlap
            )
            print(f"{k}: {goodOverlap}/{badOverlap}", end = " ")
        print()


a = []

for i in range(len(plotLevels)):
    g = [0] * 20
    b = [0] * 20
    for j in range(len(goodDash[levels.index(plotLevels[i])])):
        for k in range(20):
            g[k] += goodDash[levels.index(plotLevels[i])][j][k]
        for k in range(20):
            b[k] += badDash[levels.index(plotLevels[i])][j][k]
    g = [k/len(goodDash[levels.index(plotLevels[i])]) for k in g]
    b = [k/len(goodDash[levels.index(plotLevels[i])]) for k in b]
    entry = []
    entry += [100*sum(g[:5])/alignLenRead[levels.index(plotLevels[i])], 100*sum(b[:5])/alignLenFake[levels.index(plotLevels[i])]] 
    entry += [1000*sum(g[9:16])/alignLenRead[levels.index(plotLevels[i])], 1000*sum(b[9:16])/alignLenFake[levels.index(plotLevels[i])]]
    a.append(entry)

a = np.array(a)

fig, axs = plt.subplots(2, 1)

print(a.T)
print(a.T.shape)

axs[0].imshow(a.T[:2])
axs[0].set_xticks(np.arange(len(plotLevels)))
axs[0].set_xticklabels(plotLevels)
axs[0].set_yticks(np.arange(2))
axs[0].set_yticklabels(["<= 5 pos", "<= 5 neg"])

axs[1].imshow(a.T[2:])
axs[1].set_xticks(np.arange(len(plotLevels)))
axs[1].set_xticklabels(plotLevels)
axs[1].set_yticks(np.arange(2))
axs[1].set_yticklabels(["<= 10 pos", "<= 10 neg"])

for i in range(len(plotLevels)):
    for j in range(2):
        text = axs[0].text(i, j, str(a[i, j])[:7],
                       ha="center", va="center", color="w")
    for j in range(2, 4):
        text = axs[1].text(i, j-2, str(a[i, j])[:7],
                       ha="center", va="center", color="w")

axs[0].set_ylabel("Gap len")
axs[0].set_xlabel("Level number")
axs[1].set_ylabel("Gap len")
axs[1].set_xlabel("Level number")
fig.tight_layout()
plt.show()


'''
a = []

for i in range(len(plotLevels)):
    g = [0] * 20
    b = [0] * 20
    for j in range(len(goodDash[levels.index(plotLevels[i])])):
        for k in range(20):
            g[k] += goodDash[levels.index(plotLevels[i])][j][k]
        for k in range(20):
            b[k] += badDash[levels.index(plotLevels[i])][j][k]
    g = [k/len(goodDash[levels.index(plotLevels[i])]) for k in g]
    b = [k/len(goodDash[levels.index(plotLevels[i])]) for k in b]
    entry = [sum(g[:5])/(sum(b[:5])+0.0000001)] + [sum(g[9:16])/(sum(b[9:16])+0.0000001)]
    a.append(entry)

a = np.array(a)

fig, ax = plt.subplots()

ax.imshow(a.T)
ax.set_xticks(np.arange(len(plotLevels)))
ax.set_xticklabels(plotLevels)
ax.set_yticks(np.arange(2))
ax.set_yticklabels(["<= 5", ">= 10"])

for i in range(len(plotLevels)):
    for j in range(2):
        text = ax.text(i, j, str(a[i, j])[:7],
                       ha="center", va="center", color="w")

ax.set_ylabel("Gap len")
ax.set_xlabel("Level number")
fig.tight_layout()
plt.show()
'''


dim1, dim2 = 2, 3
fig, axs = plt.subplots(dim1, dim2)

for i in range(len(plotLevels)):
    X, Y = [], []
    for k in range(len(plotKmerLen)):
        data = overlap[levels.index(plotLevels[i])][kmerLen.index(plotKmerLen[k])]
        x, y = plotAOC(data)
        X = x
        Y.append(y)
        axs[i//dim2, i%dim2].plot(X, y, label=str(plotKmerLen[k]), linewidth=2)
    # axs[i].legend(loc="lower right")
    Y = np.array(Y)
    axs[i//dim2, i%dim2].set_title(f"Level {plotLevels[i]}")
    axs[i//dim2, i%dim2].set_aspect('equal', adjustable='box')
#axs[dim1-1, dim2-1].legend(loc="lower right")

handles, labels = axs[dim1-1, dim2-1].get_legend_handles_labels()
fig.subplots_adjust(bottom=0.1, wspace=0.1)
fig.legend(handles, labels, loc='lower center', ncol=dim1*dim2)

plt.show()

dim1, dim2 = 2, 3
fig, axs = plt.subplots(dim1, dim2)

for i in range(len(plotLevels)):
    X, Y = [], []
    for j in range(len(plotKmerLen)):
        data = pomery[levels.index(plotLevels[i])][kmerLen.index(plotKmerLen[j])]
        if len(data) == 0:
            continue
        X = range(0, 205)
        y = [0] * 205
        for k in data:
            y[math.floor(k)] += 1
        # Y.append(y)
        Y.append([sum(y[:k]) / len(data) for k in range(len(y))])
        axs[i//dim2, i%dim2].plot(list(X), Y[-1], label=str(plotKmerLen[j]), linewidth=2)
    Y = np.array(Y)
    axs[i//dim2, i%dim2].set_title(f"Level {plotLevels[i]}")

handles, labels = axs[dim1-1, dim2-1].get_legend_handles_labels()
fig.subplots_adjust(bottom=0.1, wspace=0.1)
fig.legend(handles, labels, loc='lower center', ncol=dim1*dim2)

plt.show()
