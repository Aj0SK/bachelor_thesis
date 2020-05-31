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
alignedSquiggles = "../prepareData/alignedSquiggles_1000.txt"

smoothParam = 5
repeatSignal = 10
workingLen = 2000

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
longestCommonK = [[0 for j in kmerLen] for i in levels]
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
        
        print(f"Level {levels[levels.index(l)]}: ", end = "")
        for k in kmerLen:
            counter = 0
            for i in range(len(a) - k + 1):
                w1, w2 = a[i : i + k], b[i : i + k]
                if ("-" not in w1) and ("-" not in w2):
                    counter += 1
                    break
            longestCommonK[levels.index(l)][kmerLen.index(k)] += counter
            print(f"{longestCommonK[levels.index(l)][kmerLen.index(k)]}", end = " ")
        print(f"# {readCounter}")


for l in levels:
    last_all = 0
    for i in range(len(longestCommonK[levels.index(l)])):
        if longestCommonK[levels.index(l)][i] == readNum:
            last_all = i

    print("\hline")
    print(f"{levels[levels.index(l)]}", end = " ")
    #print(f"kmerLen", end = " ")
    
    for i in range(last_all, last_all+10):
        print(f"& {kmerLen[i]}", end = " ")
    print(f"\\\\")
    #print(f"{levels[levels.index(l)]}", end = " ")
     
    for i in range(last_all, last_all+10):
        print(f"& {longestCommonK[levels.index(l)][i]}", end = " ")
    print(f"\\\\")
        
        
