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

smoothParam = 5
repeatSignal = 10
workingLen = 5000

readNum = 50

kmerLen = list(range(4, 36, 1))
levels = list(range(4, 15, 1))

plotLevels = [4, 5, 7, 9, 11, 13]  # range(4, 10)
plotKmerLen = [4, 7, 13, 17, 21, 28]  # range(4, 40, 10)

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

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

pomery = [[[] for j in kmerLen] for i in levels]
overlap = [[[] for j in kmerLen] for i in levels]
goodDash = [[] for i in levels]
badDash = [[] for i in levels]
readCounter = 0

for posRead in posReadsPaths:
    if readCounter == readNum:
        break

    nadavca_align = nadavca.align_signal(refFilePath, [posRead])

    if len(nadavca_align) != 1:
        continue

    nadavca_align = nadavca_align[0]
    fromSignal, toSignal = nadavca_align[0].signal_range

    if (toSignal - fromSignal) < workingLen:
        continue

    #print(f"Signal alligned from {fromSignal} to {toSignal}")
    print("Working on", posRead)
    print(f"So far done {readCounter} reads")
    readCounter += 1

    refSeq = "".join(nadavca_align[0].reference_part)
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


#fig, axs = plt.subplots(len(plotLevels))

import math

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
    #a.append(g)
    #a.append(b)
    #a.append([g[i]/b[i] if b[i] != 0 else 1 for i in range(len(g))])
    a.append([sum(g[:5])/(sum(b[:5])+0.0000001)] + [sum(g[9:16])/(sum(b[9:16])+0.0000001)])

plt.imshow(a, cmap = "hot", interpolation = "nearest")
plt.show()


dim1, dim2 = 2, 3
fig, axs = plt.subplots(dim1, dim2)

for i in range(len(plotLevels)):
    X, Y = [], []
    for k in range(len(plotKmerLen)):
        data = overlap[levels.index(plotLevels[i])][kmerLen.index(plotKmerLen[k])]
        x, y = plotAOC(data)
        X = x
        Y.append(y)
        axs[i//dim2, i%dim2].plot(X, y, label=str(k), linewidth=2)
    # axs[i].legend(loc="lower right")
    Y = np.array(Y)
    axs[i//dim2, i%dim2].set_title(f"This is level {plotLevels[i]}")
    axs[i//dim2, i%dim2].set_aspect('equal', adjustable='box')
axs[dim1-1, dim2-1].legend(loc="lower right")

plt.show()


dim1, dim2 = 2, 3
fig, axs = plt.subplots(dim1, dim2)

import math

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
    Y = np.array(Y)
    axs[i//dim2, i%dim2].plot(X, Y.T)
    axs[i//dim2, i%dim2].set_title(f"This is level {plotLevels[i]}")

plt.show()
