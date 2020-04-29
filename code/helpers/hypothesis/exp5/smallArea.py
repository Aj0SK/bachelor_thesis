refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
#readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

targetContig = 'contig1'
targetBeg, targetEnd = 10000, 20000

posTestCases, negTestCases = 34, 34
levels = 11
repeatSignal = 10
kmerLength = 21
overflow = 0.3
smoothParam = 5
refWindowSize = 1000
refWindowJump = 700
fromRead, toRead = 5000, 20000
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

import matplotlib
import matplotlib.pyplot as plt

def overlap(dict1, dict2):
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
    return intersect

def plotAOC(src):
    src.sort()
    src.reverse()
    X, Y = [], []
    x, y = 0, 0
    for i in src:
        X.append(x)
        Y.append(y)
        if i[1] == 1:
            y += 1
        else:
            x += 1
    print(X)
    print(Y)
    plt.scatter(X, Y)
    plt.show()

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

def f(b, e, c, d):
    if e <= c or b >= d:
        return False
    return True

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# load filenames of all positive and negative reads
posFast5 = [i for i in open("goodReads.txt", "r").readlines()]
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

assert len(posFast5) >= posTestCases, "Not enough positive testcases!"
assert len(negFast5) >= negTestCases, "Not enough negative testcases!"

posFast5 = posFast5[:posTestCases]
negFast5 = negFast5[:negTestCases]

################################################################################

hashTable = {}

for contig in Fasta(refFilePath):
    if contig.name != targetContig:
        continue
    ref = str(contig)
    ref = ref[targetBeg:targetEnd]
    contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
    hashTable = getDictFromSequence(contigSignal, refWindowSize, refWindowJump)

def processRead(path):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]

    readDict = getDictFromSequence(readSignal, refWindowSize, refWindowJump)
    hits = overlap(readDict, hashTable)
    print("Number of hits is {0}".format(hits))
    return hits

good, bad = 0, 0
helper = []
print("Positive:")

for filePath in posFast5:
    filePath = filePath.strip()
    #try:
    #    readSeq, basecallTable = getSeqfromRead(filePath)
    #except:
    #    continue
    #if len(readSeq) < (toRead // repeatSignal):
    #    continue

    #hits = [
    #    aln for aln in referenceIdx.map(readSeq)
    #    if (aln.q_en - aln.q_st > 0.95 *
    #        len(readSeq)) and aln.strand == 1 and aln.ctg==targetContig
    #    and f(aln.r_st, aln.r_en, targetBeg, targetEnd)
    #]

    #print(f"{len(hits)}")
    #if len(hits) == 0:
    #    print("Zle je")
    #    continue

    helper.append((processRead(filePath), 1))

print("\n\nNegative:")

for filePath in negFast5:
    #print(filePath)
    '''try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead // repeatSignal):
        continue
    '''
    helper.append((processRead(filePath), 0))

print("{0}/{1}".format(good, good + bad))
plotAOC(helper)
