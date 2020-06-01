refFilePath = "../../../data/sapIngB1.fa"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

readsPosFilePath = "../../../data/pos-15"
readsNegFilePath = "../../../data/neg-basecalled"

import sys

posTestCases, negTestCases = 200, 200
levels = int(sys.argv[1]) # 4
repeatSignal = 10
kmerLength = int(sys.argv[2]) # 24
overflow = 0.3
smoothParam = 5
refWindowSize = 5000
refWindowJump = 3000
fromRead, toRead = 5000, 7000
workingContig = "contig3"

################################################################################

import glob
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

import matplotlib
import matplotlib.pyplot as plt

sys.path.append("../")
from signalHelper import stringToSignal, getLevels, getSignalFromRead, getSeqfromRead, produceRandom
from signalHelper import getLevelString, computeNorm, computeString, smoothSignal, buildDictionary, overlappingKmers

def overlap(dict1, dict2):
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
            #intersect += 1*(1/dict1[kmer])
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

def plotAOC(src):
    src.sort()
    src.reverse()
    X, Y = [0], [0]
    x, y = 0, 0
    for i in src:
        if i[1] == 0:
            y += 1
        else:
            x += 1
        X.append(x)
        Y.append(y)
    return X, Y

hashTable = {}
processed = []

for contig in Fasta(refFilePath):
    if contig.name != workingContig:
        continue
    processed.append(contig.name)
    
    contigStr = str(contig)
    contigSignal = stringToSignal(contigStr, mod, repeatSignal=repeatSignal)
    hashTable = getDictFromSequence(contigSignal, refWindowSize, refWindowJump)

#print("Hashtable readyfor {0} nums!".format(contigNum))
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

def processRead(path, contigName):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]

    readString = getLevelString(readSignal, smoothParam, levels, overflow)
    readDict = buildDictionary(readString, kmerLength)
    return overlap(readDict, hashTable)

########################################

data = []

for filePath in posFast5:
    if posTestCases == 0:
        break

    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead // repeatSignal):
        continue
    hits = [
        aln for aln in referenceIdx.map(readSeq)
        if aln.q_en - aln.q_st > 0.95 * len(readSeq) and
        aln.strand == 1 and aln.ctg == workingContig
    ]
    if len(hits) == 0:
        continue
    hit = hits[0]
    
    posTestCases -= 1
    lvlStringHits = processRead(filePath, workingContig)
    data.append((lvlStringHits, 1))
    
for filePath in negFast5:
    if negTestCases == 0:
        break

    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead // repeatSignal):
        continue
    hits = [
        aln for aln in referenceIdx.map(readSeq)
        if aln.q_en - aln.q_st > 0.95 * len(readSeq) and
        aln.strand == 1 and aln.ctg == workingContig
    ]
    if len(hits) != 0:
        continue

    negTestCases -= 1
    lvlStringHits = processRead(filePath, workingContig)
    data.append((lvlStringHits, 0))
    
data.sort()

print(f"Levels: {levels} kmerLength: {kmerLength} #", end="")
for i in data:
    print(f" {i[0]} {i[1]}", end = "")
print()
