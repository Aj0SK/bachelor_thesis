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

readNum = 10

kmerLen = list(range(4, 35, 3))
levels = list(range(4, 12, 2))

plotLevels = levels  # range(4, 10)
plotKmerLen = kmerLen  # range(4, 40, 10)

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
)


def plotAOC(src):
    src.sort()
    src.reverse()
    X, Y = [], []
    x, y = 0, 0
    for i in src:
        X.append(x)
        Y.append(y)
        if i[1] == 0:
            y += 1
        else:
            x += 1
    return X, Y


################################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=0)
negReadsPaths = getReadsInFolder(readsNegFilePath, minSize=0)

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

nonZeroHits = [[0] * len(kmerLen) for _ in range(len(levels))]

pomery = [[[] for j in kmerLen] for i in levels]
overlap = [[[] for j in kmerLen] for i in levels]
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

    print(f"Signal alligned from {fromSignal} to {toSignal}")
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
        for k in kmerLen:
            counter = 0
            for i in range(len(readStrings[l]) - k + 1):
                for j in range(len(refStrings[l]) - k + 1):
                    w1 = readStrings[l][i : i + k]
                    w2 = refStrings[l][j : j + k]
                    if ("-" not in w1) and ("-" not in w2):
                        counter += 1

    for l in levels:
        # print(f"levels({l}): #", end="")
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
        print()


fig, axs = plt.subplots(len(plotLevels))

import math

for i in range(len(plotLevels)):
    X, Y = [], []
    for j in range(len(plotKmerLen)):
        if len(pomery[i][j]) == 0:
            continue
        data = pomery[i][j]
        X = range(0, 205)
        y = [0] * 205
        for k in data:
            y[math.floor(k)] += 1
        # Y.append(y)
        Y.append([sum(y[:k]) / len(data) for k in range(len(y))])
    Y = np.array(Y)
    axs[i].plot(X, Y.T)
    axs[i].set_title(f"This is level {levels[i]}")

plt.show()

fig, axs = plt.subplots(len(plotLevels))

for i in range(len(plotLevels)):
    X, Y = [], []
    for k in range(len(plotKmerLen)):
        data = overlap[levels.index(plotLevels[i])][kmerLen.index(plotKmerLen[k])]
        x, y = plotAOC(data)
        X = x
        Y.append(y)
        axs[i].plot(X, y, label=str(k))
    # axs[i].legend(loc="lower right")
    Y = np.array(Y)
    axs[i].set_title(f"This is level {levels[i]}")

plt.show()
