refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

targetContig = "contig1"
targetBeg, targetEnd = 0, 1000000000#1000000#50000

posTestCases, negTestCases = 40, 40
levels = 6
repeatSignal = 10
kmerLength = 12
overflow = 0.30
smoothParam = 5

maxCount = 5000

################################################################################

import sys
import copy
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import stringToSignal

from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
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

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

hashTable = {}

for contig in Fasta(refFilePath):
    ref = str(contig)
    contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
    contigSignal = contigSignal[targetBeg:targetEnd]
    levelStr = getLevelStr(contigSignal)
    
    for i in range(len(levelStr)-kmerLength+1):
        kmer = levelStr[i:i+kmerLength]
        hashTable[kmer] = hashTable.get(kmer, 0) + 1
    break

pocty = [0] * maxCount

for k, v in hashTable.items():
    if v >= maxCount:
        print(f"Really high count: {v} -> {k}")
        continue
    pocty[v] += 1

plt.scatter(range(1, maxCount+1), pocty, marker=r"+")
plt.axhline(y = 0)
plt.show()
