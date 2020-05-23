################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

readsPosFilePath = "../../../data/pos-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"
refFilePath = "../../../data/sapIngB1.fa"
refIndex = "../../../data/ref.index"

repeatSignal = 10
overflow = 0.3
smoothParam = 5

kmerLen = 17
level = 12

workingContig = "contig3"

readSignalBeg, readSignalEnd = 10000, 15000

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
    getSeqfromRead,
    getReadsInFolder,
    stringAllignment,
    overlappingKmers,
    buildDictionary,
    computeStringInformative,
    getLevelString,
    stringToSignal,
)
from signalHelper import computeNorm, computeString, smoothSignal

import matplotlib
import matplotlib.pyplot as plt


def overlappingKmersSpecial(str1, str2, k):
    dic1 = buildDictionary(str1, k)
    index, count = 0, 0
    while index < len(str2) - k + 1:
        w = str2[index : index + k]
        if w in dic1:
            count += 1
            index += 5
            continue
        index += 1
    return count


def buildDictionarySpecial(d, string, k):
    for i in range(0, len(string) - k + 1):
        kmer = string[i : i + k]
        if kmer not in d:
            d[kmer] = []
        d[kmer].append(i)
    return


def delWithCount(hashTable, cutoff):
    toDel = []
    for key, count in hashTable.items():
        if len(count) >= cutoff:
            toDel.append(key)
    for key in toDel:
        del hashTable[key]

def helperF(x, limit = 7):
    out = []
    x.sort()

    count, endd = 0, 0
    for begg in range(len(x) - 1):
        while endd != len(x):
            if x[begg][0] + (readSignalEnd - readSignalBeg) + 500 < x[endd][0]:
                break
            endd += 1

        if (endd - begg) <= limit:
            continue

        n = endd - begg
        dp = [0] * n
        for i in range(n):
            dp[i] = 1
            tempL = [dp[j] for j in range(0, i) if x[begg + j][1] <= x[begg + i][1]]
            if len(tempL) != 0:
                dp[i] += max(tempL)

        if max(dp) <= limit:
            continue
        
        #print(f"{x[begg:endd]}")
        
        count += 1
        out.append(x[begg][0])
    print(f"Count is {count}")
    return out

################################################################################
posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=1000000)

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)
################################################################################

globalNorms = {}
storeContig = {}
infoString = None

with open(refIndex, "r") as outFile:
    for line in outFile:
        line = line.split()
        if int(line[1]) != level or line[2] != "+":
            continue
        levelStr = line[3]
        storeContig[line[0]] = levelStr

'''for contig in ref:
    if contig.name != workingContig:
        continue
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    contigSignal = smoothSignal(contigSignal, smoothParam)
    globalNorms[contig.name] = computeNorm(contigSignal, 0, len(contigSignal))
    storeContig[contig.name], infoString = computeStringInformative(
        contigSignal,
        0,
        len(contigSignal),
        globalNorms[contig.name][0],
        globalNorms[contig.name][0],
        level,
        overflow=overflow,
    )'''

print("Preparation done!")

hashTable = {}
hashTables = []
contigNames = []

for contigName in storeContig.keys():
    hashTables.append({})
    contigNames.append(contigName)
    buildDictionarySpecial(hashTables[-1], storeContig[contigName], kmerLen)
    buildDictionarySpecial(hashTable, storeContig[contigName], kmerLen)

for sample in posReadsPaths[:400]:
    try:
        readFastq, readEvents = getSeqfromRead(sample)
    except:
        continue

    hits = [
        aln
        for aln in referenceIdx.map(readFastq)
        if aln.q_en - aln.q_st > 0.95 * len(readFastq)
        and aln.strand == 1
        and aln.ctg == workingContig
    ]
    if len(hits) != 1:
        continue
    hit = hits[0]

    location = hit.r_st / len(Fasta(refFilePath)[hit.ctg])
    print(f"I am in ctg {hit.ctg} in around {location}")

    refSeq = str(Fasta(refFilePath)[hit.ctg][hit.r_st : hit.r_en])
    refSignal = stringToSignal(refSeq, mod, repeatSignal)
    refSignal = smoothSignal(refSignal, smoothParam)
    refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
    #refShift, refScale = globalNorms[hit.ctg][0], globalNorms[hit.ctg][1]
    refString = computeString(
        refSignal, 0, len(refSignal), refShift, refScale, level, overflow=overflow
    )
    refString = refString[10 : len(refString) - 10]

    contigX = storeContig[hit.ctg]
    startInRef = -1
    for e in range(len(contigX) - len(refString) + 1):
        w = contigX[e : e + len(refString)]
        if w == refString:
            # print("Dobre")
            startInRef = e
            print(f"Hit is in ctg {hit.ctg} in [{e};{e+readSignalEnd//repeatSignal}] around {e/len(contigX)}")
    if startInRef == -1:
        print("Problem!")
        #continue

    readSignal = getSignalFromRead(sample)
    readSignal = readSignal[readSignalBeg:readSignalEnd]
    readString = getLevelString(readSignal, smoothParam, level, overflow)

    #refString = refString[readSignalBeg//repeatSignal:readSignalEnd//repeatSignal]

    u = []
    for i in range(len(readString) - kmerLen + 1):
        for j in range(len(refString) - kmerLen + 1):
            w1 = readString[i:i+kmerLen]
            w2 = refString[j:j+kmerLen]
            if w1 == w2:
                u.append((i,j))
    
    #print(u)
    #exit(0)

    refHash = {}
    for i in range(len(refString) - kmerLen + 1):
        w = refString[i : i + kmerLen]
        refHash[w] = refHash.get(w, 0) + 1

    countHits1, countHits2 = 0, 0
    hits1, hits2 = [], []
    for i in range(len(readString) - kmerLen + 1):
        w = readString[i : i + kmerLen]
        countHits1 += len(hashTable.get(w, []))
        countHits2 += min(1, refHash.get(w, 0))

        if w in hashTable:
            hits1.append(w)
        if w in refHash:
            hits2.append(w)

    #print(f"Roznych je tam {len(set(hits2))}")

    wind = 100
    x = []
    for h in range(len(hashTables)):
        if contigNames[h] != workingContig:
            continue
        contigLen = len(storeContig[contigNames[h]])
        aproxPos = int(wind * location * contigLen // contigLen)
        d = wind * [0]
        for i in range(len(readString) - kmerLen + 1):
            w = readString[i : i + kmerLen]
            hits = hashTables[h].get(w, [])
            for j in range(len(hits)):
                d[wind * hits[j] // contigLen] += 1
                x.append((hits[j], i))
        d[aproxPos] += 1000
        print(d)

        possibleHits = helperF(x, 5)
        for i in possibleHits:
            if abs(i) <= 1500 and i >= startInRef:
                print(startInRef, i)

    print(
        f"Level: {level} kmerLen: {kmerLen} hits -> {countHits1} and hits with ref only {countHits2}"
    )
