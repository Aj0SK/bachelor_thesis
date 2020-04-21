import sys
import numpy as np
import mappy as mp
from nadavca.dtw import KmerModel
from pyfaidx import Fasta
import random

sys.path.append("../")
from signalHelper import getSeqfromRead, seqSignalCor, getSignalFromRead, stringToSignal, getReadsInFolder, stringAllignment
from signalHelper import computeNorm, computeString, smoothSignal, buildDictionary, overlappingKmers
from signalHelper import countDashes

import matplotlib.pyplot as plt

##################################

refFile = "../../../data/sapIngB1.fa"
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

maxTests = 100

levels = 12
repeatSignal = 10

kmerLen = 9

signalFrom = 5000#int(sys.argv[3])
signalTo = 30000#int(sys.argv[4])

mod = KmerModel.load_from_hdf5(kmerModelFilePath)
posReads = getReadsInFolder(readsPosFilePath, minSize = 0)
negReads = getReadsInFolder(readsNegFilePath, minSize = 1000000)

totalG = 0
totalF = 0
G, F = [], []
badClas = 0

spaces_G = [0] * 14
spaces_F = [0] * 14

successfulReads = 0

rat = []

### get corresponding part of the reference using minimap2
referenceIdx = mp.Aligner(refFile)
assert referenceIdx, "failed to load/build reference index"

for readFile in posReads[:min(len(posReads), maxTests)]:
    print(readFile)
    ### read read
    try:
        readFastq, readEvents = getSeqfromRead(readFile)
    except:
        print("Bad read!")
        continue
    readSeq = seqSignalCor(signalFrom, signalTo, readEvents)
    readSignal = np.array(getSignalFromRead(readFile)[signalFrom:signalTo],
                             dtype=float)
    # readSeq - sequence cut out from read
    # readSignal - corresponding signal section

    # require a single hit with at least 95% coverage of length
    hits = [aln for aln in referenceIdx.map(readSeq)
            if aln.q_en - aln.q_st > 0.95*len(readSeq)]
    if len(hits) != 1:
        print("Too many or too few hits, skipping read.")
        continue
    hit = hits[0]
    successfulReads += 1

    if (hit.strand == 1):
        refSeq=str(Fasta(refFile)[hit.ctg][hit.r_st:hit.r_en])
        fakeSeq=str(-Fasta(refFile)[hit.ctg][hit.r_st:hit.r_en])
    else:
        refSeq=str(-Fasta(refFile)[hit.ctg][hit.r_st:hit.r_en])
        fakeSeq=str(Fasta(refFile)[hit.ctg][hit.r_st:hit.r_en])

    refSignal = np.array(stringToSignal(refSeq, mod, repeatSignal = repeatSignal),
                    float)
    #fakeSignal = np.array(stringToSignal(fakeSeq, mod, repeatSignal = repeatSignal),
    #                float)
    fakeSignal = []
    fakeIndex = -1
    while len(fakeSignal) <= signalTo:
        fakeIndex = random.randint(0, len(negReads)-1)
        fakeSignal = np.array(getSignalFromRead(negReads[fakeIndex]), dtype=float)
    fakeSignal = fakeSignal[signalFrom:signalTo]
    _, readEvents = getSeqfromRead(negReads[fakeIndex])
    fakeSeq = seqSignalCor(signalFrom, signalTo, readEvents)
    
    #print("A:", refSeq, "\n")
    #print("B:", fakeSignalSeq, "\n")

    readSignalSm = smoothSignal(readSignal,5)
    refSignalSm = smoothSignal(refSignal,5)
    fakeSignalSm = smoothSignal(fakeSignal,5)

    readShiftSm, readScaleSm = computeNorm(readSignalSm,0,len(readSignalSm))
    refShiftSm, refScaleSm = computeNorm(refSignalSm,0,len(refSignalSm))
    fakeShiftSm, fakeScaleSm = computeNorm(fakeSignalSm,0,len(fakeSignalSm))

    readString2Sm = computeString(readSignalSm,0,len(readSignalSm),readShiftSm,readScaleSm,levels,overflow=0.25)
    refString2Sm = computeString(refSignalSm,0,len(refSignalSm),refShiftSm,refScaleSm,levels,overflow=0.25)
    fakeString2Sm = computeString(fakeSignalSm,0,len(fakeSignalSm),fakeShiftSm,fakeScaleSm,levels,overflow=0.25)

    #print("readsm-1:",readString2Sm)
    #print("refsm-1 :",refString2Sm)
    #print("fakesm-1:",fakeString2Sm)
    
    print("Readstring allignment")
    a, b = stringAllignment(refString2Sm, readString2Sm)
    print(a)
    print(b)

    print("Fakestring allignment")
    c, d = stringAllignment(refString2Sm, fakeString2Sm)
    print(c)
    print(d)
    
    print("Dashes in readstring allignment")
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]:
        spaces_G[i-1] += countDashes(a, i)+countDashes(b, i)
        print(i, ":", countDashes(a, i)+countDashes(b, i))
    
    print("Fakestring in readstring allignment")
    
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]:
        spaces_F[i-1] += countDashes(c, i)+countDashes(d, i)
        print(i, ":", countDashes(c, i)+countDashes(d, i))

print("Skipped {0} out of {1}".format(maxTests-successfulReads, maxTests))

print("Good is: " + str(spaces_G))
print("Bad  is: " + str(spaces_F))
