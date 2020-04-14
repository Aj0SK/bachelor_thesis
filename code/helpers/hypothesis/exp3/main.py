# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

posTestCases, negTestCases = 1000, 10
levels = 6
repeatSignal = 10
kmerLength = 21
overflow = 0.3
smoothParam = 5
refWindowSize = 1000
refWindowJump = 700
hashTablesNum = 20
fromRead, toRead = 5000, 50000

################################################################################

import sys
import glob
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import stringToSignal, getLevels, getSignalFromRead, getSeqfromRead, produceRandom
from signalHelper import computeNorm, computeString, smoothSignal, buildDictionary, overlappingKmers

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# load filenames of all positive and negative reads
posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

assert len(posFast5) >= posTestCases, "Not enough positive testcases!"
assert len(negFast5) >= negTestCases, "Not enough negative testcases!"

################################################################################

hashTables = []
processed = []
hashWinSize = -1
hashWinJump = -1

for contig in Fasta(refFilePath):
    processed.append(contig.name)
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    hashWinSize = len(contigSignal) // hashTablesNum
    hashWinJump = int(0.95*hashWinSize)
    for i in range(0, len(contigSignal) - hashWinSize + 1, hashWinJump):
        hashTables.append({})
        for winBeg in range(i, i+hashWinSize-refWindowSize+1, refWindowJump):
            winEnd = winBeg + refWindowSize
            currSignal = copy.deepcopy(contigSignal[winBeg:winEnd])
            currSignal = np.array(currSignal, float)
            currSignal = smoothSignal(currSignal, smoothParam)
            currSignalShift, currSignalScale = computeNorm(currSignal, 0, refWindowSize)
            currString = computeString(currSignal, 0, refWindowSize, currSignalShift, currSignalScale, levels, overflow=overflow)
            currDict = buildDictionary(currString, kmerLength)
            hashTables[-1].update(currDict)
    break # only process one contig for now

print("Hashtable ready!")
#######################################
'''
print("Overlap is:")
for i in range(len(hashTables)):
    for j in range(len(hashTables)):
        if i != j:
            counter = 0
            for k in hashTables[i]:
                if k in hashTables[j]:
                    counter += 1
            print("{0} {1}: {2} with sizes {3} {4}".format(i, j, counter, len(hashTables[i]), len(hashTables[j])))
'''


def overlap(dict1, dict2):
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
    return intersect

good, bad = 0, 0

def processRead(path, goodTable=-1):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]
    readSignal = smoothSignal(readSignal, smoothParam)
    
    readDict = {}
    for winBeg in range(0, len(readSignal)-refWindowSize+1, refWindowJump):
        winEnd = winBeg + refWindowSize
        currSignal = copy.deepcopy(readSignal[winBeg:winEnd])
        currSignal = np.array(currSignal, float)
        currSignalShift, currSignalScale = computeNorm(currSignal, 0,
                                                   len(currSignal))
        currString = computeString(currSignal,
                               0,
                               refWindowSize,
                               currSignalShift,
                               currSignalScale,
                               levels,
                               overflow=overflow)

        readDict.update(buildDictionary(currString, kmerLength))
    hits = [overlap(readDict, hashTable) for hashTable in hashTables]
    max_hits = max(hits)
    max_i = -1
    
    for i in range(len(hits)):
        if i == goodTable:
            print("->", end='')
        if hits[i] == max_hits:
            print("max->", end='')
            max_i = i
        print(hits[i], end=' ')
    print()
    global good, bad
    if goodTable != -1 and max_i == goodTable:
        print("Good one!")
        good += 1
    if goodTable != -1 and max_i != goodTable:
        bad += 1
    return


########################################

print("Positive:")

for filePath in posFast5[:posTestCases]:
    #print(filePath)
    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead//repeatSignal):
        continue
    hits = [
        aln for aln in referenceIdx.map(readSeq)
        if (aln.q_en - aln.q_st > 0.95 * len(readSeq))
        and aln.strand == 1 and aln.ctg in processed
    ]
    if len(hits) != 1:
        continue
    table_num = (repeatSignal * hits[0].r_st) // hashWinJump
    processRead(filePath, table_num)

print("\n\nNegative:")

for filePath in negFast5[:negTestCases]:
    #print(filePath)
    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead//repeatSignal):
        continue
    processRead(filePath)

print("{0}/{1}".format(good, good+bad))
