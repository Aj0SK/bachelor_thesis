refFilePath = "../../../data/sapIngB1.fa"
kmerModelFilePath = "../../../data/kmer_model.hdf5"
refIndex = "../../../data/ref.index"

levels = 4
kmerLength = 14
overflow = 0.30
smoothParam = 5

maxCount = 1000000

################################################################################

import sys
import copy
import numpy as np
from pyfaidx import Fasta
import matplotlib
import matplotlib.pyplot as plt

################################################################################

hashTable = {}

with open(refIndex, 'r') as outFile:
    for line in outFile:
        line = line.split()
        if int(line[1]) != levels:
            continue
        
        levelStr = line[3]
        
        for i in range(len(levelStr)-kmerLength+1):
            kmer = levelStr[i:i+kmerLength]
            hashTable[kmer] = hashTable.get(kmer, 0) + 1

pocty = [0] * maxCount
totalNum = 0

for k, v in hashTable.items():
    totalNum += v
    if v >= maxCount:
        print(f"Really high count: {v} -> {k}")
        continue
    pocty[v-1] += 1

sumUpTo = []
suma = 0
for i in range(maxCount):
    suma += pocty[i] * (i+1)
    if abs(suma/totalNum-1.0) < 0.01:
        break
    sumUpTo.append(suma/totalNum)

plt.xlabel("Početnosť kmerov")
plt.ylabel("Agregované pokrytie referencie")
plt.plot(list(range(1, len(sumUpTo)+1)), sumUpTo)
plt.axhline(y = 0)
plt.axhline(y = 1)
plt.show()
