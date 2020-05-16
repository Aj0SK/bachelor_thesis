refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

repeatSignal = 10
overflow = 0.30
smoothParam = 5
fromRead, toRead = 10000, 20000
kmerNum = 300000
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

def buildDictionarySpecial(d, string, k):
    for i in range(0, len(string) - k + 1):
        kmer = string[i : i + k]
        if kmer not in d:
            d[kmer] = 0
        d[kmer] += 1
    return

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

posReads = getReadsInFolder(readsPosFilePath, minSize = 0)
negReads = getReadsInFolder(readsNegFilePath, minSize = 0)

################################################################################

def helper(hashTable, reads, infoString):
    counter, total, totalNum = 0, 0, 0

    for readFile in reads:
        if totalNum == kmerNum:
            break
        try:
            readFastq, _ = getSeqfromRead(readFile)
        except:
            continue

        hits = [
            aln
            for aln in referenceIdx.map(readFastq)
            if aln.q_en - aln.q_st > 0.95 * len(readFastq) and 
            aln.strand == 1
        ]

        if infoString == "+"  and len(hits) != 1:
            continue
        if infoString == "-"  and len(hits) != 0:
            continue

        readSignal = getSignalFromRead(readFile)
            
        if len(readSignal) <= toRead:
            continue
            
        counter += 1
            
        readSignal = readSignal[fromRead:toRead]
        readString = getLevelStr(readSignal, levels)
                
        before = total
        for i in range(len(readString)):
            if kmerNum == totalNum:
                break
            kmer = readString[i : i + k]
            totalNum += 1
            total += hashTable.get(kmer, 0)
        
        print(f"Read {readFile} : \t {total-before}")
    print(f"{infoString} k {k} l {levels} -> {total} / {totalNum}")

for levels in [8]:#range(6, 15):
    storeContig = {}

    for contig in Fasta(refFilePath):
        print("Next!")
        ref = str(contig)
        contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
        levelStr = getLevelStr(contigSignal, levels)
        storeContig[contig.name] = levelStr

    #print("Refstrings ready!")

    hk = list(range(5, 35))

    for k in hk:
        hashTable = {}
        
        for contig in storeContig.values():
            buildDictionarySpecial(hashTable, contig, k)
        #toDel = []
        #for key, count in hashTable.items():
        #    if count >= 100000:
        #        toDel.append(key)
        #for key in toDel:
        #    del hashTable[key]
        
        helper(hashTable, posReads, "+")
        helper(hashTable, negReads, "-")
