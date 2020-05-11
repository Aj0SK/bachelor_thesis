refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

targetBeg, targetEnd = 0, 10000000

posTestCases, negTestCases = 40, 40
repeatSignal = 10
overflow = 0.30
smoothParam = 5
fromRead, toRead = 5000, 20000
kmerNum = 10000
################################################################################

import sys
import random
import glob
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import (
    stringToSignal,
    getLevels,
    getSignalFromRead,
    getSeqfromRead,
    getReadsInFolder,
    produceRandom,
)
from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
    overlappingKmers,
)

import matplotlib
import matplotlib.pyplot as plt


def getLevelStr(signal, l):
    currSignal = np.array(copy.deepcopy(signal), float)
    currSignal = smoothSignal(currSignal, smoothParam)
    currSignalShift, currSignalScale = computeNorm(currSignal, 0, len(currSignal))
    currString = computeString(
        currSignal,
        0,
        len(currSignal),
        currSignalShift,
        currSignalScale,
        l,
        overflow=overflow,
    )
    return currString

def buildDictionary(dict, string, k):
    for i in range(0, len(string) - k + 1):
        kmer = string[i : i + k]
        if kmer not in dict:
            dict[kmer] = 0
        dict[kmer] += 1
    return

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

posReads = getReadsInFolder(readsPosFilePath)
negReads = getReadsInFolder(readsNegFilePath)

################################################################################

for levels in range(4, 15):
    storeContig = {}

    for contig in Fasta(refFilePath):
        ref = str(contig)[:targetEnd]
        contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
        levelStr = getLevelStr(contigSignal, levels)
        storeContig[contig.name] = levelStr

    print("Refstrings ready!")

    hk = list(range(5, 35))

    for k in hk:

        hashTable = {}

        for contig in storeContig.values():
            buildDictionary(hashTable, contig, k)
        
        toDel = []
        for key, count in hashTable.items():
            if count >= 10000:
                toDel.append(key)
        
        for key in toDel:
            del hashTable[key]
        
        total, counter, totalNum = 0, 0, 0

        for readFile in posReads:
            if counter == posTestCases or totalNum == kmerNum:
                break
            try:
                readFastq, _ = getSeqfromRead(readFile)
            except:
                continue

            hits = [
                aln
                for aln in referenceIdx.map(readFastq)
                if aln.q_en - aln.q_st > 0.95 * len(readFastq) and aln.strand == 1
            ]

            if len(hits) != 1:
                continue
                
                #print("Bingo!")

            readSignal = getSignalFromRead(readFile)[fromRead:toRead]
            readString = getLevelStr(readSignal, levels)
                
            counter += 1
                
            for i in range(len(readString)):
                if kmerNum == totalNum:
                    break
                kmer = readString[i : i + k]
                totalNum += 1
                total += hashTable.get(kmer, 0)
            #print("Counter +1!")
            #print(f"Add {len(readString)}")
        
        print(f"+ k {k} l {levels} -> {total} / {totalNum}")
        totalNum, total, counter = 0, 0, 0
        for readFile in negReads:
            if counter == negTestCases or totalNum == kmerNum:
                break
            try:
                readFastq, _ = getSeqfromRead(readFile)
            except:
                continue

            hits = [
                aln
                for aln in referenceIdx.map(readFastq)
                if aln.q_en - aln.q_st > 0.10 * len(readFastq) and aln.strand == 1
            ]

            if len(hits) != 0:
                continue
            
            readSignal = getSignalFromRead(readFile)[fromRead:toRead]
            readString = getLevelStr(readSignal, levels)
                
            counter += 1
                
            for i in range(len(readString)):
                if kmerNum == totalNum:
                    break;
                kmer = readString[i : i + k]
                totalNum += 1
                total += hashTable.get(kmer, 0)
        print(f"- k {k} l {levels} -> {total} / {totalNum}")
