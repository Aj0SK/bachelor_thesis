refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "goodReadsDebug.txt"
readsNegFilePath = "../../../data/neg-basecalled"

targetContig = "contig1"
targetBeg, targetEnd = 0, 1000000#50000

posTestCases, negTestCases = 40, 40
levels = 6
repeatSignal = 10
kmerLength = 12
overflow = 0.30
smoothParam = 5
refWindowSize = 1000
refWindowJump = 700
fromRead, toRead = 5000, 20000
contigNum = 1

################################################################################

import sys
import random
import glob
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import (
    stringToSignal,
    getLevels,
    getSignalFromRead,
    getSeqfromRead,
    produceRandom,
)
from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
    buildDictionary,
    overlappingKmers,
)

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


def buildDictionarySpecial(string, k):
    dict = {}
    for i in range(0, len(string) - k + 1):
        kmer = string[i : i + k]
        if kmer not in dict:
            dict[kmer] = []
        dict[kmer].append(i)
    return dict


def getDictFromSequence(signal, l = False):
    dic = {}
    currSignal = np.array(copy.deepcopy(signal), float)
    currSignal = smoothSignal(currSignal, smoothParam)
    currSignalShift, currSignalScale = computeNorm(currSignal, 0, len(currSignal))
    currString = computeString(
        currSignal,
        0,
        len(currSignal),
        currSignalShift,
        currSignalScale,
        levels,
        overflow=overflow,
    )
    if l == True:
        return currString
    return buildDictionarySpecial(currString, kmerLength)


def intervalOverlap(b, e, c, d):
    if e <= c or b >= d:
        return False
    return True


################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# load filenames of all positive and negative reads
data = [i.strip() for i in open(readsPosFilePath, "r").readlines() if i.strip() != ""]
posFast5 = [data[i] for i in range(0, len(data), 2)]
basecalledFast5 = [data[i] for i in range(1, len(data), 2)]

negFast5 = glob.glob(readsNegFilePath + "/*.fast5", recursive=True)

################################################################################

hashTable = {}

for contig in Fasta(refFilePath):
    if contig.name != targetContig:
        continue
    ref = str(contig)
    ref = ref[targetBeg:targetEnd]
    contigSignal = stringToSignal(ref, mod, repeatSignal=repeatSignal)
    hashTable = getDictFromSequence(contigSignal)
    break

#for k in sorted(hashTable, key=hashTable.get, reverse=True)[:100]:
#    pass
#    # print("{0} {1}".format(k, hashTable[k]))
#    # del hashTable[k]


def processRead(path, readFromRef=False):
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]
    
    readString = getDictFromSequence(readSignal, l = True)
    
    #hits = overlap(readDict, hashTable)
    #return hits

    myHits = 0
    jump = 2*kmerLength
    operateRange = 10*kmerLength
    for i in range(100):
        beg = random.randint(0, len(readString)-operateRange+1)
        w = readString[beg:beg+kmerLength]
        candidates = hashTable.get(w, [])

        while len(candidates) > 1:
            j = random.randint(beg+jump, beg+operateRange-kmerLength+1)

            w = readString[j:j+kmerLength]
            newcand = []
            for k in candidates:
                for l in hashTable.get(w, []):
                    if abs(k-l) < operateRange*1.30:
                        newcand.append(k)
                        break

            candidates = newcand


        if len(candidates) != 0:
            myHits += 1
            print(f"Candidates len is {len(candidates)}")
            print(str(candidates)[:200])

    print(f"My hits is {myHits}")
    #print("Number of hits is {0}".format(hits))

    return myHits


helper = []
print("Positive:")

for i in range(len(posFast5)):
    filePath = posFast5[i]
    readSeq = basecalledFast5[i]
    hits = [
        aln
        for aln in referenceIdx.map(readSeq)
        if (aln.q_en - aln.q_st > 0.95 * len(readSeq))
        and aln.strand == 1
        and aln.ctg == targetContig
        and intervalOverlap(aln.r_st, aln.r_en, targetBeg, targetEnd)
    ]

    if len(hits) == 0:
        negTestCases -= 1
        continue

    if not (hits[0].r_st > targetBeg and hits[0].r_en < targetEnd):
        negTestCases -= 1
        continue
    
    print(100*"#")
    print(hits[0].r_st)
    helper.append((processRead(filePath, readFromRef=True), 1))

print("\n\nNegative:")

exit(0)

for i in range(negTestCases):
    print(100*"#")
    filePath = negFast5[i]
    helper.append((processRead(filePath), 0))

plotAOC(helper)
