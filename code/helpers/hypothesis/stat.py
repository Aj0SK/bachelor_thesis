# reference file that we use to generate signal
refFilePath = "../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../data/pos-basecalled"
readsNegFilePath = "../../data/neg-basecalled"

maxTests = 100

# what part of signal we take from read
fromSignal, toSignal = 11000, 13000
refWindowSize = 100
refWindowJump = 1

# number of levels we use
levels = 9
repeatSignal = 10
smoothParam = 5
overflow = 0.1
kmerLength = 5

import sys
import math
import random
import h5py
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import (
    getSeqfromRead,
    seqSignalCor,
    getSignalFromRead,
    stringToSignal,
    getReadsInFolder,
)
from signalHelper import (
    computeNorm,
    computeString,
    smoothSignal,
    buildDictionary,
    overlappingKmers,
)

################################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

posReads = getReadsInFolder(readsPosFilePath)
negReads = getReadsInFolder(readsNegFilePath)

################################################################################


merStat = {}

def getDictFromSequence(signal, refWindowSize, refWindowJump):
    dic = {}
    for winBeg in range(0, len(signal) - refWindowSize + 1, refWindowJump):
        winEnd = winBeg + refWindowSize
        currSignal = np.array(copy.deepcopy(signal[winBeg:winEnd]), float)
        currSignal = smoothSignal(currSignal, smoothParam)
        currSignalShift, currSignalScale = computeNorm(currSignal, 0, refWindowSize)
        currString = computeString(
            currSignal,
            0,
            refWindowSize,
            currSignalShift,
            currSignalScale,
            levels,
            overflow=overflow,
        )
        newDict = buildDictionary(currString, kmerLength)
        for i in newDict:
            dic[i] = dic.get(i, 0) + newDict[i]
    return dic


for readFile in posReads[:maxTests]:
    try:
        readFastq, readEvents = getSeqfromRead(readFile)
    except:
        print("Bad read!")
        continue
    # require a single hit with at least 95% coverage of length
    hits = [
        aln
        for aln in referenceIdx.map(readFastq)
        if aln.q_en - aln.q_st > 0.95 * len(readFastq) and aln.strand == 1
    ]
    if len(hits) != 1:
        #print("Too many or too few hits, skipping read.")
        continue
    hit = hits[0]

    refStr = str(ref[hit.ctg][hit.r_st : hit.r_en])

    refSignal = stringToSignal(refStr, mod, repeatSignal=repeatSignal)
    readSignal = getSignalFromRead(readFile)[fromSignal:toSignal]
    fakeSignal = []
    while len(fakeSignal) <= toSignal:
        fakeSignal = getSignalFromRead(negReads[random.randint(0, len(negReads) - 1)])
    fakeSignal = fakeSignal[fromSignal:toSignal]

    refSignal = np.array(refSignal, float)
    readSignal = np.array(readSignal, float)
    fakeSignal = np.array(fakeSignal, float)

    readSignal = smoothSignal(readSignal,5)
    refSignal = smoothSignal(refSignal,5)
    fakeSignal = smoothSignal(fakeSignal,5)
    
    refShift, refScale = computeNorm(refSignal, 0, len(refSignal))
    readShift, readScale = computeNorm(readSignal, 0, len(readSignal))
    fakeShift, fakeScale = computeNorm(fakeSignal, 0, len(fakeSignal))
    
    readString = computeString(readSignal, 0, len(readSignal), readShift, readScale , levels, overflow = 0.3)
    refString = computeString(refSignal, 0, len(refSignal), refShift, refScale, levels, overflow = 0.3)
    fakeString = computeString(fakeSignal, 0, len(fakeSignal), fakeShift, fakeScale, levels, overflow = 0.3)

    for i in range(0, len(readString)-3+1):
        mer = readString[i:i+3]
        merStat[mer] = merStat.get(mer, 0) + 1

for k in sorted(merStat, key=merStat.get, reverse=True)[:30]:
    print("{0} {1}".format(k, merStat[k]))
