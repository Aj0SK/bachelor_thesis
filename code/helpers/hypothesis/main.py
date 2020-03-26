# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../data/pos-basecalled"
readsNegFilePath = "../../data/neg-basecalled"

testCases = 150
levels = 5
kmerLength = 17
refWindowSize = 1000
refJump = 300

################################################################################

import sys
import copy
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

from signalHelper import stringToSignal, getLevels, getSignalFromRead
from signalHelper import computeNorm, computeString, smoothSignal, buildDictionary, overlappingKmers


################################################################################
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

hashTable = {}

for contig in Fasta(refFilePath):
    print("Next contig!")
    contigSignal = stringToSignal(str(contig), mod, repeatSignal = 10)
    print(len(contigSignal))
    for winBeg in range(0, len(contigSignal)-refWindowSize+1, refJump):
        winEnd = winBeg + refWindowSize
        currSignal = copy.deepcopy(contigSignal[winBeg:winEnd])
        currSignal = np.array(currSignal, float)
        currSignal = smoothSignal(currSignal, 5)
        currSignalShift, currSignalScale = computeNorm(currSignal, 0, refWindowSize)
        currString = computeString(currSignal, 0, refWindowSize, currSignalShift, currSignalScale, levels)
        currDict = buildDictionary(currString, kmerLength)
        hashTable.update(currDict)
    #break    
    #out = getLevels(np.array(contigSignal, dtype = float))

    #for i in out:
    #    hashTable[i] = 1
print("Hashtable ready!")
print(str(hashTable)[:1000])
########################################
### hashTable helper

# for a whole read in form of levelStrings return total number of hits in *hTable*
def countMatch(hTable, signalStrings):
    counter = 0
    for i in signalStrings:
        if i in hTable:
            counter += 1
    return counter


def processRead(path):
    #print("\nProcessing " + path)
    readSignal = np.array(getSignalFromRead(path), dtype = float)
    if readSignal.shape[0] < 9000:
        return
    
    readSignal = readSignal[8000:9000]
    readSignal = smoothSignal(readSignal, 5)
    readSignalShift, readSignalScale = computeNorm(readSignal, 0, len(readSignal))
    readString = computeString(readSignal, 0, len(readSignal), readSignalShift, readSignalScale, levels)
    hits = countMatch(hashTable, buildDictionary(readString, kmerLength))
    print(hits, end = ' ')
    return

########################################

import glob

# load filenames of all positive and negative reads
posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

assert len(posFast5) >= testCases, "Not enough positive testcases!"
assert len(negFast5) >= testCases, "Not enough negative testcases!"

print("Positive:")

for filePath in posFast5[:testCases]:
    processRead(filePath)

print()
print("Negative:")

for filePath in negFast5[:testCases]:
    processRead(filePath)

print()
