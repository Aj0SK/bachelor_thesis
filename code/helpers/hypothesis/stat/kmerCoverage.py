refFilePath = "../../../data/sapIngB1.fa"
refIndex = "../../../data/ref.index"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

targetBeg, targetEnd = 0, 1000000000#1000000#50000

kmerLens = [4, 7, 13, 17, 21, 28] # list(range(4, 36, 1))
levels = [4, 5, 7, 9, 11, 13] # list(range(4, 15, 1))
maxCount = 1000000

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

def helper(hashTable):
    pocty = [0] * (maxCount+1)
    totalNum = 0

    for k, v in hashTable.items():
        totalNum += v
        if v >= maxCount:
            continue
        pocty[v] += 1
        
    sumUpTo = []
    suma = 0
    for i in range(maxCount):
        suma += pocty[i] * (i+1)
        if abs(suma/totalNum-1.0) < 0.01:
            break
        sumUpTo.append(suma/totalNum)

    return sumUpTo

################################################################################

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

hashTable = {}

storeContig = {}
contigs = []

with open(refIndex, "r") as outFile:
    for line in outFile:
        line = [l.strip() for l in line.split()]
        if int(line[1]) not in levels or line[2] != "+":
            continue

        levelStr = line[3]
        storeContig[(int(line[1]), line[0])] = levelStr
        contigs.append(line[0])

contigs = set(contigs)
hashTable = {}
results = {}

for l in levels:
    for k in kmerLens:
        h = {}
        for contig in contigs:
            levelStr = storeContig[(l, contig)]
            for i in range(len(levelStr)-k+1):
                kmer = levelStr[i:i+k]
                h[kmer] = h.get(kmer, 0) + 1
            #hashTable[(l, k)] = h
        results[(l, k)] = helper(h)
        print(f"l: {l} k: {k}")

dim1, dim2 = 2, 3
fig, axs = plt.subplots(dim1, dim2)

for i in range(len(levels)):
    axs[i//dim2, i%dim2].axhline(y = 0)
    axs[i//dim2, i%dim2].axhline(y = 1)
    for k in kmerLens:
        ee = results[(levels[i], k)]
        #plt.xlabel("Početnosť kmerov")
        #plt.ylabel("Agregované pokrytie referencie")
        axs[i//dim2, i%dim2].plot(list(range(1, len(ee)+1)), ee)
plt.show()
