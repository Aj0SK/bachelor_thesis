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

kmerLens = [34] #list(range(31, 35, 1)) # list(range(30, 39, 1)) # list(range(4, 20, 1))
levels = [4] #[12]  # [4, 5, 6, 7, 8, 9, 10, 11, 12]

workingLen = 10000

readSignalBeg, readSignalEnd = 7000, 15000

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
from signalHelper import computeNorm, computeStringInformative, computeString, smoothSignal

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

globalNorms = {}


globalNorms = {'contig1': (-0.07255694315479783, 1.1418830906171593), 'contig2': (-0.07693958064866857, 1.142595822123931), 'contig3': (-0.0733159470348583, 1.1403456974048618), 'contig4': (-0.07498020221643004, 1.1460200046328326), 'contig5': (-0.07389364528397992, 1.145021580548445), 'mtDNA': (0.012250303979524157, 1.144555008498658)}
'''for contig in ref:
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    contigSignal = smoothSignal(contigSignal, smoothParam)
    globalNorms[contig.name] = computeNorm(contigSignal, 0, len(contigSignal))

print(str(globalNorms))'''

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
        
        hashTable = {}
        hashTables = {}
        contigNames = []

        for contigName in storeContig.keys():
            hashTables[contigName] = {}
            contigNames.append(contigName)
            #buildDictionarySpecial(hashTables[contigName], storeContig[contigName], kmerLen)
            buildDictionarySpecial(hashTable, storeContig[contigName], kmerLen)

        print("Reads:")
        for sample in posReadsPaths[:100]:

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

            refPosition = hit.r_st/len(ref[hit.ctg])
            print(f"I am in ctg {hit.ctg} in around {refPosition}")

            if hit.strand == 1:
                refSeq = str(ref[hit.ctg][hit.r_st : hit.r_en])

            refSignal = stringToSignal(refSeq, mod, repeatSignal)
            refSignal = smoothSignal(refSignal, smoothParam)
            #refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
            refShift, refScale = globalNorms[hit.ctg][0], globalNorms[hit.ctg][1]
            refString = computeString(
                refSignal,
                0,
                len(refSignal),
                refShift,
                refScale,
                level,
                overflow=overflow,
            )
            refString = refString[2:-2]

            readSignal = getSignalFromRead(sample)
            readSignalLen = len(readSignal)
            readSignal = readSignal[readSignalBeg:readSignalEnd]
            readString = getLevelString(readSignal, smoothParam, level, overflow)

            found = None
            print(f"Len of refstring is {len(refString)}")
            for i in range(len(storeContig[hit.ctg])-len(refString) + 1):
                w = storeContig[hit.ctg][i:i+len(refString)]
                if w == refString:
                    found = i
                    print(f"Nasiel som ho v {i}")

            if found == None:
                print("Problem")
                exit(0)

            refString = refString[:int(len(refString)*((readSignalEnd+2000)/readSignalLen))]

            refHash = {}
            for i in range(len(refString) - kmerLen + 1):
                w = refString[i : i + kmerLen]
                refHash[w] = refHash.get(w, 0) + 1

            countHits1, countHits2 = 0, 0
            hits1, hits2 = [], []
            for i in range(len(readString) - kmerLen + 1):
                w = readString[i : i + kmerLen]

                if w in hashTable:
                    countHits1 += len(hashTable.get(w, []))
                    hits1.append(w)
                if w in refHash:
                    countHits2 += min(1, refHash.get(w, 0))
                    hits2.append(w)
            
            '''good = False
            for w in hits2:
                for index in hashTables[hit.ctg].get(w, []):
                        if good == True:
                            break
                        if abs(index-found) < len(refString):
                            print(f"Juchuchu {(index/len(storeContig[hit.ctg]))}")
                            good = True
                            break'''

            overlap = sum([1 for w in hits1 if w in hits2])

            print(
                f"Level: {level} kmerLen: {kmerLen} hits -> {countHits1} and hits with ref only {countHits2} with {overlap}"
            )
