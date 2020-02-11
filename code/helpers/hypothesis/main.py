refFilePath = "../../data/sapIngB1.fa"
kmerModelFilePath = "../../data/kmer_model.hdf5"

################################################################################

import sys
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

from signalHelper import stringToSignal, getLevels

################################################################################
genes = Fasta(refFilePath)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

genesSignal = {}

for i in genes:
    genesSignal[i.name] = stringToSignal(str(i), mod)

hashTable = {}

for key in genesSignal.keys():
    out = getLevels(np.array(genesSignal[key], dtype = float))
    for i in out:
        if i not in hashTable:
            hashTable[i] = 0
        hashTable[i] += 1

print("Hashtable ready!")
########################################
### hashTable

import h5py
import glob

def getRead(filename):
    readFile = h5py.File(filename, 'r')
    readName = str(list(readFile['Raw/Reads'])[0])
    rawData = readFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    return rawData

def countMatch(hTable, signalStrings):
    counter = 0
    for i in signalStrings:
        if i in hTable:
            counter += 1
    return counter

########################################

readsPosFilePath = "../../data/pos"
readsNegFilePath = "../../data/neg"

posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

print("Positive:")

for filePath in posFast5:
    posRead = np.array(getRead(filePath), dtype = float)
    oneReadLevels = getLevels(posRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(posRead)), end=' ')

print()
print("Negative:")

for filePath in negFast5:
    negRead = np.array(getRead(filePath), dtype = float)
    oneReadLevels = getLevels(negRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(negRead)), end=' ')

print()
