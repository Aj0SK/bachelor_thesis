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

maxTests = 400
repeatSignal = 10

kmerLen = list(range(4, 40))
levels = list(range(4, 10))

plotLevels = range(4, 10)
plotKmerLen = range(4, 40, 10)

signalFrom = 5000
signalTo = 7000


def helper(a, b):
    a = list(a)
    b = list(b)

    for s in range(1, len(a) - 1):
        if b[s - 1] != "-" and b[s + 1] != "-" and a[s + 1] != "-" and a[s - 1] != "-":
            if a[s] == "-" and ord(b[s]) == (ord(a[s - 1]) + ord(a[s + 1])) // 2:
                a[s] = "X"  # chr((ord(a[s-1]) + ord(a[s+1]))//2)
            if b[s] == "-" and ord(a[s]) == (ord(b[s - 1]) + ord(b[s + 1])) // 2:
                b[s] = "X"  # chr((ord(a[s-1]) + ord(a[s+1]))//2)

    a = "".join(a)
    b = "".join(b)
    return a, b


def succesive(a, b):
    longest = 0
    currLen = 0

    for i in range(len(a)):
        if a[i] == "-" or b[i] == "-":
            longest = max(longest, currLen)
            currLen = 0
        else:
            currLen += 1
    return max(longest, currLen)


rat = []

### get corresponding part of the reference using minimap2
referenceIdx = mp.Aligner(refFile)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)
posReads = getReadsInFolder(readsPosFilePath, minSize=1000000)
negReads = getReadsInFolder(readsNegFilePath, minSize=1000000)

goodK, totalK = 0, 0

dobre, zle = 0, 0

pomery = [[[] for j in kmerLen] for i in levels]

successfulReads = 0

for readFile in posReads[: min(len(posReads), maxTests)]:
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
    successfulReads += 1

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
    

print("Skipped {0} out of {1}".format(maxTests - successfulReads, maxTests))

print("Good is: " + str(spaces_G))
print("Bad  is: " + str(spaces_F))
