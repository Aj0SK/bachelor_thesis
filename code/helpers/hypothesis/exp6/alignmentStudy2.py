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

repeatSignal = 10

kmerLen = list(range(4, 40))
levels = list(range(4, 10))

plotLevels = range(4, 10)
plotKmerLen = range(4, 40, 10)

workingLen = 2000

minSignal, maxSignal = -2.0, 2.0

import sys
import math
import random
import h5py
import numpy as np
import mappy as mp
from pyfaidx import Fasta
import nadavca
from nadavca.dtw import KmerModel

sys.path.append("../")

from signalHelper import (
    stringToSignal,
    getSignalFromRead,
    getReadsInFolder,
    stringAllignment,
    overlappingKmers,
    buildDictionary,
)
from signalHelper import computeNorm, computeString, smoothSignal

import matplotlib
import matplotlib.pyplot as plt

def overlappingKmersSpecial(str1, str2, k):
    dic1 = buildDictionary(str1, k)
    index, count = 0, 0
    while index < len(str2)-k+1:
        w = str2[index:index+k]
        if w in dic1:
            count += 1
            index += 5
            continue
        index += 1
    return count
        

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
    #print(X)
    #print(Y)
    return X, Y

############s###################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=0)
negReadsPaths = getReadsInFolder(readsNegFilePath, minSize=0)

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

nonZeroHits = [[0] * len(kmerLen) for _ in range(len(levels))]

pomery = [[[] for j in kmerLen] for i in levels]
readCounter = 0

for posRead in posReadsPaths:
    nadavca_align = nadavca.align_signal(
        refFilePath, [posRead]
    )

    if len(nadavca_align) != 1:
        continue
    if nadavca_align[0][0].reverse_complement == True:
        #print("Problem!")
        pass
    
    if readCounter == 150:
        break
    
    nadavca_align = nadavca_align[0]

    # refStr = str(ref[nadavca_align[0].contig_name])
    fromSignal, toSignal = nadavca_align[0].signal_range
    
    if (toSignal-fromSignal) < workingLen:
        continue
    
    print("Working on", posRead)
    print(f"So far done {readCounter} reads")
    readCounter += 1
    
    refSeq = "".join(nadavca_align[0].reference_part)

    refSignal = np.array(stringToSignal(refSeq, mod, repeatSignal=repeatSignal), float)
    readSignal = np.array(getSignalFromRead(posRead)[fromSignal:toSignal], dtype=float)
    # fakeSignal = np.array(stringToSignal(fakeSeq, mod, repeatSignal = repeatSignal),
    #                float)
    fakeSignal = []
    fakeIndex = -1
    while len(fakeSignal) <= toSignal:
        fakeIndex = random.randint(0, len(negReadsPaths) - 1)
        fakeSignal = np.array(getSignalFromRead(negReadsPaths[fakeIndex]), dtype=float)
    fakeSignal = fakeSignal[fromSignal:toSignal]

    readSignal = readSignal[:workingLen]
    refSignal = refSignal[:workingLen]
    fakeSignal = fakeSignal[:workingLen]

    readSignal = smoothSignal(readSignal, 5)
    refSignal = smoothSignal(refSignal, 5)
    fakeSignal = smoothSignal(fakeSignal, 5)
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
    
    for l in levels[:1]:
        a, b = stringAllignment(refStrings[l], readStrings[l])
        c, d = stringAllignment(refStrings[l], fakeStrings[l])
        a, b, c, d = a[0:120], b[0:120], c[0:120], d[0:120]
        #print(f"{l}:\n{a}\n{b}\n{c}\n{d}\n")
        
    '''for l in levels:
        for k in kmerLen:
            goodOverlap = overlappingKmers(refStrings[l], readStrings[l], k)
            badOverlap = overlappingKmers(refStrings[l], fakeStrings[l], k)
            print(f"{l}-levels, {k}-kmerlen: {goodOverlap} vs {badOverlap}")'''
            
    
    for l in levels:
        print(f"levels({l}): #", end = "")
        for k in kmerLen:
            #goodOverlap = overlappingKmersSpecial(refStrings[l], readStrings[l], k)
            #badOverlap = overlappingKmersSpecial(refStrings[l], fakeStrings[l], k)
            goodOverlap = overlappingKmers(refStrings[l], readStrings[l], k)
            badOverlap = overlappingKmers(refStrings[l], fakeStrings[l], k)
            #if goodOverlap >= 3:
            #    nonZeroHits[l-levels[0]][k-kmerLen[0]] += 1
            #x = nonZeroHits[l-levels[0]][k-kmerLen[0]]
            #if k >= l and k < 30:
                #print(f"{k}:{x}", end = '\t')
            print(f"{k}:{goodOverlap}/{badOverlap}", end = '\t')
            
            pomery[levels.index(l)][kmerLen.index(k)].append((goodOverlap, 0))
            pomery[levels.index(l)][kmerLen.index(k)].append((badOverlap, 1))
            
            if goodOverlap == 0:
                badOverlap = 1000
                goodOverlap = 100
            if (badOverlap/goodOverlap) > 2.0:
                badOverlap = 2*goodOverlap
                #rat.append(100.0 * badOverlap / goodOverlap)
            #pomery[l-levels[0]][k-kmerLen[0]].append(100.0 * badOverlap / goodOverlap)
        print()
            
'''
fig, axs = plt.subplots(len(levels))

import math

for i in range(len(levels)):https://www.forbes.com/sites/qai/2020/05/18/boeing-and-delta-among-stocks-to-watch-this-week/#44277ebe1f64
    #for j in range(len(kmerLen)):
    X, Y = [], []
    for j in range(len(kmerLen)):
        if len(pomery[i][j]) == 0:
            continue
        rat = pomery[i][j]
        X = range(0, 205)
        y = [0] * 205
        for k in rat:
            y[math.floor(k)] += 1
        Y.append([sum(y[:k])/len(rat) for k in range(len(y))])
    Y = np.array(Y)
    axs[i].plot(X, Y.T)
    axs[i].set_title(f"This is level {levels[i]}")

plt.show()
'''

fig, axs = plt.subplots(len(plotLevels))

for i in range(len(plotLevels)):
    X, Y = [], []
    for k in range(len(plotKmerLen)):
        rat = pomery[levels.index(plotLevels[i])][kmerLen.index(plotKmerLen[k])]
        x, y = plotAOC(rat)
        X = x
        Y.append(y)
        axs[i].plot(X, y, label = str(k))
    axs[i].legend(loc="lower right")
    Y = np.array(Y)

plt.show()
