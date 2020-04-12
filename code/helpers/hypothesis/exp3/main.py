# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

posTestCases, negTestCases = 1000, 400
levels = 14
repeatSignal = 10
kmerLength = 25
refWindowSize = 1000
refWindowJump = 700
hashTablesNum = 5
fromRead, toRead = 8000, 14000

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
    hashWinJump = int(1.0*hashWinSize)
    for i in range(0, len(contigSignal) - hashWinSize + 1, hashWinJump):
        print("Juchuchu1 " + str(i))
        ''' for winBeg in range(i, i+hashLen-refWindowSize+1, refJump):
            winEnd = winBeg + refWindowSize
            currSignal = copy.deepcopy(contigSignal[winBeg:winEnd])
            currSignal = np.array(currSignal, float)
            currSignal = smoothSignal(currSignal, 5)
            currSignalShift, currSignalScale = computeNorm(currSignal, 0, refWindowSize)
            currString = computeString(currSignal, 0, refWindowSize, currSignalShift, currSignalScale, levels, overflow=0.3)
            currDict = buildDictionary(currString, kmerLength)
            hashTables[i//hashLen].update(currDict)'''
        currSignal = copy.deepcopy(contigSignal[i:i + hashWinSize])
        currSignal = np.array(currSignal, float)
        currSignal = smoothSignal(currSignal, 5)
        currSignalShift, currSignalScale = computeNorm(currSignal, 0,
                                                       hashWinSize)
        currString = computeString(currSignal,
                                   0,
                                   hashWinSize,
                                   currSignalShift,
                                   currSignalScale,
                                   levels,
                                   overflow=0.3)
        hashTables.append(buildDictionary(currString, kmerLength))
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

def processRead(path, goodTable=-1):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    #if readSignal.shape[0] < toRead:
    #    return
    #readSignal = readSignal[fromRead:toRead]
    readSignal = smoothSignal(readSignal, 5)
    readSignalShift, readSignalScale = computeNorm(readSignal, 0,
                                                   len(readSignal))
    readString = computeString(readSignal,
                               0,
                               len(readSignal),
                               readSignalShift,
                               readSignalScale,
                               levels,
                               overflow=0.3)

    readDict = buildDictionary(readString, kmerLength)
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
    if max_i == goodTable:
        print("Good one!")
    return


########################################

print("Positive:")

for filePath in posFast5[:posTestCases]:
    #print(filePath)
    readSeq, basecallTable = getSeqfromRead(filePath)
    hits = [
        aln for aln in referenceIdx.map(readSeq)
        if (aln.q_en - aln.q_st > 0.95 * len(readSeq))
        and aln.strand == 1 and aln.ctg in processed
    ]
    if len(hits) != 1:
        continue
    table_num = (repeatSignal * hits[0].r_st) // hashWinJump
    print("In table {0}".format(table_num))
    processRead(filePath, table_num)

print("\n\nNegative:")

for filePath in negFast5[:negTestCases]:
    #print(filePath)
    processRead(filePath)
