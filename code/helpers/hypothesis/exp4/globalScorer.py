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
    buildDictionary,
    overlappingKmers,
)
from signalHelper import countDashes

import matplotlib.pyplot as plt

refFile = "../../../data/sapIngB1.fa"
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

maxTests = 400

kmerLen = 17
levels = 14
repeatSignal = 10

signalFrom = 5000
signalTo = 7000

totalG = 0
totalF = 0
G, F = [], []
badClas = 0

spaces_G = [0] * 14
spaces_F = [0] * 14

successfulReads = 0


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

for readFile in posReads[: min(len(posReads), maxTests)]:
    #print(readFile)
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
    print(readFile)
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
    readString2Sm = computeString(
        readSignalSm,
        0,
        len(readSignalSm),
        readShiftSm,
        readScaleSm,
        levels,
        overflow=0.25,
    )
    refString2Sm = computeString(
        refSignalSm, 0, len(refSignalSm), refShiftSm, refScaleSm, levels, overflow=0.25
    )
    fakeString2Sm = computeString(
        fakeSignalSm,
        0,
        len(fakeSignalSm),
        fakeShiftSm,
        fakeScaleSm,
        levels,
        overflow=0.25,
    )

    # print("readsm-1:",readString2Sm)
    # print("refsm-1 :",refString2Sm)
    # print("fakesm-1:",fakeString2Sm)

    over1 = overlappingKmers(readString2Sm, refString2Sm, kmerLen)
    over2 = overlappingKmers(fakeString2Sm, refString2Sm, kmerLen)
    print("Overlap of len {0} -> good:{1} vs fake:{2}".format(kmerLen, over1, over2))

    a, b = stringAllignment(refString2Sm, readString2Sm)
    c, d = stringAllignment(refString2Sm, fakeString2Sm)

    oldGoodK = goodK
    for i in range(0, len(a) - kmerLen + 1):
        totalK += 1
        goodK += 1
        for j in range(i, i + kmerLen):
            if a[j] == "-" or b[j] == "-":
                goodK -= 1
                break

    if goodK != oldGoodK:
        dobre += 1
        print(f"Pocet {kmerLen}-ic je {goodK-oldGoodK}")
    else:
        zle += 1

    print(f"Pocet zarovnani so zhodujucou {kmerLen}-ticou: {dobre} vs zle: {zle}")

    # print(f"Len of readstring is {len(readString2Sm)}")
    # print(f"Len of ref is {len(refString2Sm)}")
    # print(f"Len of allignment is {len(a)}")

    a, b = helper(a, b)
    c, d = helper(c, d)

    s1 = succesive(a, b)
    s2 = succesive(c, d)

    print(f"Longest in ref vs fake: {s1} vs {s2}")

    # print("Dashes in readstring allignment")

    spaces1 = [countDashes(a[30:], i) + countDashes(b[30:], i) for i in range(1, 15)]
    spaces2 = [countDashes(c[30:], i) + countDashes(d[30:], i) for i in range(1, 15)]

    print("Number of continuos spaces in allignment of len >= 10: ", spaces1[9:])
    print("Number of continuos spaces in allignment of len >= 10: ", spaces2[9:])

    # print(f"Bigger that 10: {sum(spaces1[9:])}")
    # print(f"Bigger that 10: {sum(spaces2[9:])}")

    for i in range(1, 15):
        spaces_G[i - 1] += spaces1[i - 1]
        #print(i, ":", spaces1[i - 1])

    # print("Fakestring in readstring allignment")

    for i in range(1, 15):
        spaces_F[i - 1] += spaces2[i - 1]
        #print(i, ":", spaces2[i - 1])

    """
    Print alligned strings
    fromAlli, toAlli = 0, 200
    a, b, c, d = (
        a[fromAlli:toAlli],
        b[fromAlli:toAlli],
        c[fromAlli:toAlli],
        d[fromAlli:toAlli],
    )
    print("Readstring allignment")
    print(a)
    print(b)
    print("Fakestring allignment")
    print(c)
    print(d)
    """
    
    # print(f"Toto nas zaujima -> {goodK}/{totalK}")

print("Skipped {0} out of {1}".format(maxTests - successfulReads, maxTests))

print("Good is: " + str(spaces_G))
print("Bad  is: " + str(spaces_F))
