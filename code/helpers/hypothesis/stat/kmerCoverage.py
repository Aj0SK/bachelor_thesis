refFilePath = "../../../data/sapIngB1.fa"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

targetBeg, targetEnd = 0, 1000000000#1000000#50000

levels = 9
repeatSignal = 10
kmerLength = 14
overflow = 0.30
smoothParam = 5

maxCount = 10000

################################################################################

import sys
import copy
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
import matplotlib
import matplotlib.pyplot as plt

sys.path.append("../")
from signalHelper import (
    stringToSignal,
    computeNorm,
    computeString,
    smoothSignal,
)

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

pocty = [0] * (maxCount+1)
totalNum = 0

for k, v in hashTable.items():
    totalNum += v
    if v >= maxCount:
        print(f"Really high count: {v} -> {k}")
        continue
    pocty[v] += 1

plt.scatter(range(0, maxCount+1), pocty, marker=r"+")
plt.axhline(y = 0)
plt.show()

sumUpTo = []
suma = 0
for i in range(1, maxCount+1):
    suma += pocty[i] * i
    sumUpTo.append(suma/totalNum)

plt.xlabel("Početnosť kmerov")
plt.ylabel("Agregované pokrytie referencie")
plt.scatter(range(1, maxCount+1), sumUpTo, marker=r"+")
plt.axhline(y = 0)
plt.show()
