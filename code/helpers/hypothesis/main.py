# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../data/pos"
readsNegFilePath = "../../data/neg"

################################################################################

import sys
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

from signalHelper import stringToSignal, getLevels

################################################################################
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

hashTable = {}

for contig in Fasta(refFilePath):
    contigSignal = stringToSignal(str(contig), mod)
    out = getLevels(np.array(contigSignal, dtype = float))

    for i in out:
        hashTable[i] = 1

print("Hashtable ready!")
########################################
### hashTable

import h5py
import glob

# return signal as list from single read
def getSignalFromRead(filename):
    readFile = h5py.File(filename, 'r')
    readName = str(list(readFile['Raw/Reads'])[0])
    rawData = readFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    return rawData

# for a whole read in form of levelStrings return total number of hits in *hTable*
def countMatch(hTable, signalStrings):
    counter = 0
    for i in signalStrings:
        if i in hTable:
            counter += 1
    return counter

########################################

# load filenames of all positive and negative reads
posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

print("Positive:")

for filePath in posFast5:
    posRead = np.array(getSignalFromRead(filePath), dtype = float)
    oneReadLevels = getLevels(posRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(posRead)), end=' ')

print()
print("Negative:")

for filePath in negFast5:
    negRead = np.array(getSignalFromRead(filePath), dtype = float)
    oneReadLevels = getLevels(negRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(negRead)), end=' ')

print()
