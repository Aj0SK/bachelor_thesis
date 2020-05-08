refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "goodReadsDebug.txt"
readsNegFilePath = "../../../data/neg-basecalled"

targetContig = "contig1"
targetBeg, targetEnd = 0, 1000000#50000

posTestCases, negTestCases = 40, 40
levels = 6
repeatSignal = 10
kmerLength = 12
overflow = 0.30
smoothParam = 5
refWindowSize = 1000
refWindowJump = 700
fromRead, toRead = 5000, 20000
contigNum = 1

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
    produceRandom,
)
from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
    buildDictionary,
    overlappingKmers,
)

import matplotlib
import matplotlib.pyplot as plt

def getLevelStr(signal):
    currSignal = np.array(copy.deepcopy(signal), float)
    currSignal = smoothSignal(currSignal, smoothParam)
    currSignalShift, currSignalScale = computeNorm(currSignal, 0, len(currSignal))
    currString = computeString(
        currSignal,
        0,
        len(currSignal),
        currSignalShift,
        currSignalScale,
        levels,
        overflow=overflow,
    )
    return currString

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

hashTable = {}

for contig in Fasta(refFilePath):
    ref = str(contig)
    contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
    levelStr = getLevelStr(contigSignal)
    
    for i in range(len(levelStr)-kmerLength+1):
        kmer = levelStr[i:i+kmerLength]
        hashTable[kmer] = hashTable.get(kmer, 0) + 1
    break

pocty = [0] * 40

for k, v in hashTable.items():
    if v >= len(pocty):
        print(f"Mega: {v} -> {k}")
        continue
    pocty[v] += 1
    
for i in range(len(pocty)):
    print(f"{i}: {pocty[i]}")

plt.plot(pocty)
plt.show()
