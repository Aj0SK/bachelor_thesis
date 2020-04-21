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
hashTablesNum = 10
hashWinSize = 500000
hashWinJump = int(0.95*hashWinSize)
fromRead, toRead = 5000, 50000
contigNum = 1

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

def overlap(dict1, dict2):
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
    return intersect

def getDictFromSequence(signal, refWindowSize, refWindowJump):
    dic = {}
    for winBeg in range(0, len(signal) - refWindowSize + 1, refWindowJump):
        winEnd = winBeg + refWindowSize
        currSignal = np.array(
            copy.deepcopy(signal[winBeg:winEnd]), float)
        currSignal = smoothSignal(currSignal, smoothParam)
        currSignalShift, currSignalScale = computeNorm(currSignal, 0,
                                                       refWindowSize)
        currString = computeString(currSignal,
                                   0,
                                   refWindowSize,
                                   currSignalShift,
                                   currSignalScale,
                                   levels,
                                   overflow=overflow)
        dic.update(buildDictionary(currString, kmerLength))
    return dic


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

hashTables = {}
processed = []

for contig in Fasta(refFilePath):
    processed.append(contig.name)
    hashTables[contig.name] = []
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    for i in range(0, len(contigSignal) - hashWinSize + 1, hashWinJump):
        hashTables[contig.name].append(
            getDictFromSequence(contigSignal[i:i+hashWinSize],
                                refWindowSize,
                                refWindowJump))
    contigNum -= 1
    if contigNum == 0:
        break

print("Hashtable readyfor {0} nums!".format(contigNum))
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

#######################################
good, bad = 0, 0

def processRead(path, contigName, goodTable=-1):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]

    readDict = getDictFromSequence(readSignal, refWindowSize, refWindowJump)
    hits = [
        overlap(readDict, hashTable) for hashTable in hashTables[contigName]
    ]
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
        print("Good match!")
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
    if len(readSeq) < (toRead // repeatSignal):
        continue
    hits = [
        aln for aln in referenceIdx.map(readSeq)
        if (aln.q_en - aln.q_st > 0.95 *
            len(readSeq)) and aln.strand == 1 and aln.ctg in processed
    ]
    if len(hits) != 1:
        continue
    table_num = (repeatSignal * hits[0].r_st) // hashWinJump
    processRead(filePath, hits[0].ctg, table_num)

print("\n\nNegative:")

for filePath in negFast5[:negTestCases]:
    #print(filePath)
    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead // repeatSignal):
        continue
    for conName in processed:
        print("For contig " + conName)
        processRead(filePath, conName)

print("{0}/{1}".format(good, good + bad))
