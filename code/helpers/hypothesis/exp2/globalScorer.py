import sys
import numpy as np
import mappy as mp
from nadavca.dtw import KmerModel
from pyfaidx import Fasta
import random

sys.path.append("../")
from signalHelper import getSeqfromRead, seqSignalCor, getSignalFromRead, stringToSignal, getReadsInFolder
from signalHelper import computeNorm, computeString, smoothSignal, buildDictionary, overlappingKmers

import matplotlib.pyplot as plt

##################################

refFile = "../../../data/sapIngB1.fa"
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

maxTests = 500

levels = 9
repeatSignal = 10

kmerLen = 19

signalFrom = 5000#int(sys.argv[3])
signalTo = 20000#int(sys.argv[4])

mod = KmerModel.load_from_hdf5(kmerModelFilePath)
posReads = getReadsInFolder(readsPosFilePath, minSize = 0)
negReads = getReadsInFolder(readsNegFilePath, minSize = 1000000)

totalG = 0
totalF = 0
G, F = [], []
badClas = 0

successfulReads = 0

rat = []

### get corresponding part of the reference using minimap2
referenceIdx = mp.Aligner(refFile)
assert referenceIdx, "failed to load/build reference index"

for readFile in posReads[:min(len(posReads), maxTests)]:
    print(readFile)
    ### read read
    readFastq, readEvents = getSeqfromRead(readFile)
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
    fakeSignal = np.array(stringToSignal(fakeSeq, mod, repeatSignal = repeatSignal),
                    float)
    #fakeSignal = []
    #while len(fakeSignal) <= signalTo:
    #    fakeSignal = np.array(getSignalFromRead(negReads[random.randint(0, len(negReads)-1)]), dtype=float)
    #fakeSignal = fakeSignal[signalFrom:signalTo]
    
    #print(readSeq)
    #print(refSeq)
    #print(fakeSeq)
    # refSeq - part of the reference sequence corresponding to the read segment

    readSignalSm = smoothSignal(readSignal,5)
    refSignalSm = smoothSignal(refSignal,5)
    fakeSignalSm = smoothSignal(fakeSignal,5)

    readShift, readScale = computeNorm(readSignal,0,len(readSignal))
    refShift, refScale = computeNorm(refSignal,0,len(refSignal))
    readShiftSm, readScaleSm = computeNorm(readSignalSm,0,len(readSignalSm))
    refShiftSm, refScaleSm = computeNorm(refSignalSm,0,len(refSignalSm))
    fakeShiftSm, fakeScaleSm = computeNorm(fakeSignalSm,0,len(fakeSignalSm))

    readString = computeString(readSignal,0,len(readSignal),readShift,readScale,levels)
    readString2 = computeString(readSignal,0,len(readSignal),readShift,readScale,levels,overflow=0.25)
    readStringSm = computeString(readSignalSm,0,len(readSignalSm),readShiftSm,readScaleSm,levels)
    readString2Sm = computeString(readSignalSm,0,len(readSignalSm),readShiftSm,readScaleSm,levels,overflow=0.25)
    
    refString = computeString(refSignal,0,len(refSignal),refShift,refScale,levels)
    refString2 = computeString(refSignal,0,len(refSignal),refShift,refScale,levels,overflow=0.25)
    refStringSm = computeString(refSignalSm,0,len(refSignalSm),refShiftSm,refScaleSm,levels)
    refString2Sm = computeString(refSignalSm,0,len(refSignalSm),refShiftSm,refScaleSm,levels,overflow=0.25)
    
    fakeString2Sm = computeString(fakeSignalSm,0,len(fakeSignalSm),fakeShiftSm,fakeScaleSm,levels,overflow=0.25)

    '''
    print("read-0:",readString)
    print("ref-0 :",refString)
    print("read-1:",readString2)
    print("ref-1 :",refString2)
    print("readsm-0:",readStringSm)
    print("refsm-0 :",refStringSm)
    print("readsm-1:",readString2Sm)
    print("refsm-1 :",refString2Sm)
    print("fakesm-1:",fakeString2Sm)

    print("Dlzka:",len(readString2Sm))
    print("read vs. reference")
    '''
    x = overlappingKmers(readString2Sm,refString2Sm,kmerLen)
    totalG += x
    print(kmerLen,":",x,end="\n")
    
    y = overlappingKmers(readString2Sm,fakeString2Sm,kmerLen)
    totalF += y
    print(kmerLen,":",y,end="\n")
    
    G.append(x)
    F.append(y)
    
    print(totalG)
    print(totalF)
    if x<y:
        badClas += 1
    
    if x != 0 and ((y/x)<2.0):
        rat.append(100.0*y/x)
    
    print(totalG/totalF)
    print("Unsuccessful reads to total reads: {0}/{1}".format(badClas, successfulReads))

print("Skipped {0} out of {1}".format(maxTests-successfulReads, maxTests))

rat = np.array(rat)
fig, axs = plt.subplots()

# We can set the number of bins with the `bins` kwarg
axs.hist(rat, bins=200)

plt.show()
