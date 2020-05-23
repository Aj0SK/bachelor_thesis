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
refIndex = "../../../data/ref.index"

repeatSignal = 10
overflow = 0.3
smoothParam = 5

kmerLens = [17]  # list(range(4, 20, 1))
levels = [12]  # [4, 5, 6, 7, 8, 9, 10, 11, 12]

workingLen = 10000

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
    getLevelString,
    stringToSignal,
)
from signalHelper import computeNorm, computeString, smoothSignal

import matplotlib
import matplotlib.pyplot as plt


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


################################################################################

posReadsPaths = getReadsInFolder(readsPosFilePath, minSize=1000000)
negReadsPaths = getReadsInFolder(readsNegFilePath, minSize=1000000)

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

ref = Fasta(refFilePath)

# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

# sample = posReadsPaths[7]
"""
globalNorms = {}

for contig in ref:
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    contigSignal = smoothSignal(contigSignal, smoothParam)
    globalNorms[contig.name] = computeNorm(contigSignal, 0, len(contigSignal))
"""
print("Preparation done!")

for level in levels:
    storeContig = {}

    with open(refIndex, "r") as outFile:
        for line in outFile:
            line = line.split()
            if int(line[1]) != level or line[2] != "+":
                continue

            levelStr = line[3]
            storeContig[line[0]] = levelStr

    for kmerLen in kmerLens:
        print("Reads:")
        hashTable = {}
        hashTables = []
        contigNames = []

        for contigName in storeContig.keys():
            hashTables.append({})
            contigNames.append(contigName)
            # buildDictionarySpecial(hashTables[-1], storeContig[contigName], kmerLen)
            buildDictionarySpecial(hashTable, storeContig[contigName], kmerLen)

        # delWithCount(hashTable, 1000)

        for sample in posReadsPaths[:30]:

            try:
                readFastq, readEvents = getSeqfromRead(sample)
            except:
                continue

            hits = [
                aln
                for aln in referenceIdx.map(readFastq)
                if aln.q_en - aln.q_st > 0.95 * len(readFastq) and aln.strand == 1
            ]
            if len(hits) != 1:
                # print("Fail!")
                continue
            hit = hits[0]

            print(f"It is in ctg {hit.ctg}")
            # print(f"I am in ctg {hit.ctg} in around {hit.r_en/len(Fasta(refFilePath)[hit.ctg])}")

            if hit.strand == 1:
                refSeq = str(Fasta(refFilePath)[hit.ctg][hit.r_st : hit.r_en])
            # else:
            #    refSeq = str(-Fasta(refFilePath)[hit.ctg][hit.r_st : hit.r_en])

            refSignal = stringToSignal(refSeq, mod, repeatSignal)
            refSignal = smoothSignal(refSignal, smoothParam)
            refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
            # refShift, refScale = globalNorms[hit.ctg][0], globalNorms[hit.ctg][1]
            refString = computeString(
                refSignal,
                0,
                len(refSignal),
                refShift,
                refScale,
                level,
                overflow=overflow,
            )

            readSignal = getSignalFromRead(sample)
            readSignal = readSignal[readSignalBeg:readSignalEnd]
            readString = getLevelString(readSignal, smoothParam, level, overflow)

            refHash = {}
            for i in range(len(refString) - kmerLen + 1):
                w = refString[i : i + kmerLen]
                refHash[w] = refHash.get(w, 0) + 1
                # if w not in hashTable:
                #    print("Problem")

            countHits1, countHits2 = 0, 0
            hits1, hits2 = [], []
            for i in range(len(readString) - kmerLen + 1):
                w = readString[i : i + kmerLen]

                # if len(hashTable.get(w, [])) > 100:
                #    continue

                countHits1 += len(hashTable.get(w, []))
                countHits2 += min(1, refHash.get(w, 0))

                if w in hashTable:
                    hits1.append(w)
                if w in refHash:
                    hits2.append(w)
                # hitsSpecial += refHash.get(w, 0)

            # if countHits1 > 1000000:
            #    continue

            """for h in range(len(hashTables)):
                hity = 0
                d = [0]*10
                for i in range(len(readString) - kmerLen + 1):
                    w = readString[i : i + kmerLen]
                    hits = hashTables[h].get(w, [])
                    if len(hits) > 100:
                        continue
                    hity += len(hits)
                    for pom in hits:
                        d[10*pom//len(storeContig[contigNames[h]])] += 1
                        
                print(f"{h} : {hity} -> {d}")"""

            # countHits2 = overlappingKmers(refString, readString, kmerLen)

            overlap = sum([1 for w in hits1 if w in hits2])

            print(
                f"Level: {level} kmerLen: {kmerLen} hits -> {countHits1} and hits with ref only {countHits2} with {overlap}"
            )
