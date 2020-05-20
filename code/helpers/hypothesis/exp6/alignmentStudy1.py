import sys
import numpy as np
import mappy as mp
from nadavca.dtw import KmerModel
from pyfaidx import Fasta
import random

sys.path.append("../")
from signalHelper import (
    getSeqfromRead,
    seqSignalCor,
    getSignalFromRead,
    stringToSignal,
    getReadsInFolder,
    stringAllignment,
)
from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
    overlappingKmers,
)
from signalHelper import countDashes

import matplotlib.pyplot as plt

refFile = "../../../data/sapIngB1.fa"
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

maxTests = 100
repeatSignal = 10

levels = list(range(4, 12))
kmerLen = list(range(4, 40))

plotLevels = list(range(4, 12))
plotKmerLen = (range(4, 40))

signalFrom = 5000
signalTo = 7000

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

### get corresponding part of the reference using minimap2
referenceIdx = mp.Aligner(refFile)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)
posReads = getReadsInFolder(readsPosFilePath, minSize=1000000)
negReads = getReadsInFolder(readsNegFilePath, minSize=1000000)

pomery = [[[] for j in kmerLen] for i in levels]

for readFile in posReads:
    print(readFile)
    try:
        readFastq, readEvents = getSeqfromRead(readFile)
    except:
        continue
    readSeq = seqSignalCor(signalFrom, signalTo, readEvents)

    hits = [
        aln
        for aln in referenceIdx.map(readSeq)
        if aln.q_en - aln.q_st > 0.95 * len(readSeq)
    ]
    if len(hits) != 1:
        # print("Too many or too few hits, skipping read.")
        continue
    hit = hits[0]

    maxTests -= 1
    if maxTests == 0:
        break

    if hit.strand == 1:
        refSeq = str(Fasta(refFile)[hit.ctg][hit.r_st : hit.r_en])
        fakeSeq = str(-Fasta(refFile)[hit.ctg][hit.r_st : hit.r_en])
    else:
        refSeq = str(-Fasta(refFile)[hit.ctg][hit.r_st : hit.r_en])
        fakeSeq = str(Fasta(refFile)[hit.ctg][hit.r_st : hit.r_en])

    readSignal = np.array(getSignalFromRead(readFile)[signalFrom:signalTo], dtype=float)
    refSignal = np.array(stringToSignal(refSeq, mod, repeatSignal=repeatSignal), float)
    # fakeSignal = np.array(stringToSignal(fakeSeq, mod, repeatSignal = repeatSignal),
    #                float)
    fakeSignal = []
    fakeIndex = -1
    while len(fakeSignal) <= signalTo:
        fakeIndex = random.randint(0, len(negReads) - 1)
        fakeSignal = np.array(getSignalFromRead(negReads[fakeIndex]), dtype=float)
    fakeSignal = fakeSignal[signalFrom:signalTo]

    readSignalSm = smoothSignal(readSignal, 5)
    refSignalSm = smoothSignal(refSignal, 5)
    fakeSignalSm = smoothSignal(fakeSignal, 5)
    readShiftSm, readScaleSm = computeNorm(readSignalSm, 0, len(readSignalSm))
    refShiftSm, refScaleSm = computeNorm(refSignalSm, 0, len(refSignalSm))
    fakeShiftSm, fakeScaleSm = computeNorm(fakeSignalSm, 0, len(fakeSignalSm))

    readStrings, refStrings, fakeStrings = {}, {}, {}

    for l in levels:
        for k in kmerLen:
            readStrings[l] = computeString(
                readSignalSm,
                0,
                len(readSignalSm),
                readShiftSm,
                readScaleSm,
                l,
                overflow=0.25,
            )
            refStrings[l] = computeString(
                refSignalSm,
                0,
                len(refSignalSm),
                refShiftSm,
                refScaleSm,
                l,
                overflow=0.25,
            )
            fakeStrings[l] = computeString(
                fakeSignalSm,
                0,
                len(fakeSignalSm),
                fakeShiftSm,
                fakeScaleSm,
                l,
                overflow=0.25,
            )
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
