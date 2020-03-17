# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../data/pos-basecalled"
readsNegFilePath = "../../data/neg-basecalled"

maxFakeReads = 6

# number of levels we use
numLevels = 7
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 8
kernelLen = 7
winSize = 15 * kernelLen

import math
import h5py
import copy
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import getLevels, getLevelString, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead, Table_Iterator, seqSignalCor, getReadsInFolder

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

################################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

posReads = getReadsInFolder(readsPosFilePath)
negReads = getReadsInFolder(readsNegFilePath)

################################################################################
def examineRead(sampleRead):
    fromSignal, toSignal = 11000, 60000
    fastqSeq, basecallTable = getSeqfromRead(sampleRead)
    originalSignal = np.array(getSignalFromRead(sampleRead), dtype = float)
    originalSignal = originalSignal[fromSignal:toSignal]

    signalFrTo = seqSignalCor(fromSignal, toSignal, basecallTable)
    hits = [i for i in sequenceIndex.map(signalFrTo) if i.q_en - i.q_st > 90*len(signalFrTo)//100]
    if len(hits) != 1:
        print("Too many hits")
        return
    print("Mapped from {0} to {1}".format(hits[0].q_st, hits[0].q_en))
    signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])
    artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)

    levelO = getLevels(originalSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
    levelO = set(levelO)
    levelA = getLevels(artifSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)

    # count have many strings are shared by artifSignal and originalSignal
    counterO = len([i for i in levelA if i in levelO])

    print("Artificial match rate with original signal {0}/{1}".format(counterO, len(levelO)))

    # load original signal from read
    
    fakeCounter = 0
    for j in negReads:
        fakeCounter += 1
        if fakeCounter == maxFakeReads:
            break
        fromSignalFake, toSignalFake = 11000, 60000
        fakeSignal = getSignalFromRead(j)
        fakeSignal = np.array(fakeSignal, dtype = float)
        fakeSignal = fakeSignal[fromSignalFake:toSignalFake]
        levelF = getLevels(fakeSignal, kernelLen = kernelLen, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
        counterF = len([i for i in levelA if i in levelF])
        print("Total fake match rate with original signal {0}/{1}".format(counterF, len(levelF)))

    return

################################################################################
# load basecalled sequence and signal

for i in posReads:
    print("Working on " + i)
    examineRead(i)
