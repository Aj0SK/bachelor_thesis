refFile = "../../../data/sapIngB1.fa"

# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

levels = 4
repeatSignal = 10
overflow = 0.3
smoothParam = 5
kmerLen = 19

import sys
import glob
import numpy as np
from pyfaidx import Fasta
import h5py
import mappy as mp
from nadavca.dtw import KmerModel

sys.path.append("../")
from signalHelper import getSeqfromRead, stringToSignal, getSignalFromRead, getReadsInFolder
from signalHelper import computeNorm, computeString, smoothSignal
from signalHelper import stringAllignment, countDashes
from signalHelper import overlappingKmers

def getGlobalString(signal):
    signal = smoothSignal(signal, smoothParam)
    signalShift, signalScale = computeNorm(signal, 0, len(signal))
    signalString = computeString(signal,
                                 0,
                                 len(signal),
                                 signalShift,
                                 signalScale,
                                 levels,
                                 overflow=overflow)
    return signalString


referenceIdx = mp.Aligner(refFile)
assert referenceIdx, "failed to load/build reference index"

ref = Fasta(refFile)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

posReadsFiles = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negReadsFiles = getReadsInFolder(readsNegFilePath)

contig = ref[0]
refSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
refString = getGlobalString(refSignal)
helperDict = {"a": "A", "b": "C", "c": "G", "d": "T"}
#refString = "".join(helperDict[i] for i in refString)

#with open("./nieco.fa", "w") as f:
#    f.write(">" + "x" + "\n")
#    f.write(refString.upper() + "\n")

refLevelIdx = mp.Aligner("./helper.fa")
assert refLevelIdx, "failed to load/build reference index"

print("Priprava hotovo!")

posCounter, negCounter = 3, 3
goodPosReads, counter = 0, 0

print("Positive reads!")

for readFile in posReadsFiles:
    try:
        readFastq, readEvents = getSeqfromRead(readFile)
    except:
        print("Bad read!")
        continue

    # require a single hit with at least 95% coverage of length
    hits = [
        aln for aln in referenceIdx.map(readFastq)
        if aln.q_en - aln.q_st > 0.95 *
        len(readFastq) and aln.strand == 1 and (aln.ctg == "contig1")
    ]
    if len(hits) != 1:
        #print("Too many or too few hits, skipping read.")
        continue

    counter += 1
    print(readFile)

    readSignal = np.array(getSignalFromRead(readFile), float)
    readString = getGlobalString(readSignal)

    readString = "".join(helperDict[i] for i in readString)
    readString = readString[:1500]
    
    levelHits = list(refLevelIdx.map(readString, cs = True))
    levelHits = [i for i in levelHits if i.strand == 1]
    
    if len(levelHits) == 0:
        print("No hits!")
        continue
    
    if abs((levelHits[0].r_st/len(refString))-(hits[0].r_st/len(contig)))<0.001:
        goodPosReads += 1

    hitLen = levelHits[0].r_en-levelHits[0].r_st

    print("Hit len is {0}".format(hitLen))

    print("Goodreads to all reads {0}/{1}".format(goodPosReads, counter))
    
    '''
    for i in levelHits:
        print("Hit len is {0} from {1} to {2}".format(i.r_en-i.r_st, i.q_st, i.q_en))
        print("Len of reference is {0}".format(len(contig)))
        print("Len of refstring is {0}".format(len(refString)))
        print("Hit je v {0}".format(hits[0].r_st/len(contig)))
        print("Hit je v {0}".format(i.r_st/len(refString)))
    counter += 1
    '''
    a, b = stringAllignment(refString[levelHits[0].r_st:levelHits[0].r_en], readString[levelHits[0].q_st:levelHits[0].q_en])
    
    for i in range(1, 20):
        print(i, ":", countDashes(a, i)+countDashes(b, i))
    
    
    #x = overlappingKmers(refString[levelHits[0].r_st:levelHits[0].r_en],readString[levelHits[0].q_st:levelHits[0].q_en],kmerLen)
    #print("Overlaps is {0}".format(x))
    if counter >= posCounter:
        break

print("\n\nNegative reads!\n")
print("*"*200)

for readFile in negReadsFiles[:negCounter]:
    print(readFile)
    try:
        readFastq, readEvents = getSeqfromRead(readFile)
    except:
        print("Bad read!")
        continue
    
    hits = [
        aln for aln in referenceIdx.map(readFastq)
        if aln.strand == 1 and (aln.ctg == "contig1")
    ]
    
    print("Aj pri negativnom je pocet hitov {0}".format(len(hits)))
    
    readSignal = np.array(getSignalFromRead(readFile), float)
    readString = getGlobalString(readSignal)
    readString = "".join(helperDict[i] for i in readString)
    
    print(len(readString))
    readString = readString[:1500]
    
    levelHits = list(refLevelIdx.map(readString))
    
    print("Number of hits is {0}".format(len(levelHits)))
    
    levelHits = [i for i in levelHits if i.strand == 1]
    
    if len(levelHits) == 0:
        print("Nula hitov vo fake stringu.")
        continue
    
    hitLen = levelHits[0].r_en-levelHits[0].r_st
    
    for i in levelHits:
        print("Hit len je {0}".format(i.r_en-i.r_st))
        print("Hit je v {0}".format(i.r_st/len(refString)))
    
    a, b = stringAllignment(refString[levelHits[0].r_st:levelHits[0].r_en], readString[levelHits[0].q_st:levelHits[0].q_en])
    
    #x = overlappingKmers(refString[levelHits[0].r_st:levelHits[0].r_en],readString[levelHits[0].q_st:levelHits[0].q_en],kmerLen)
    #print("Overlaps is {0}".format(x))
    
    for i in range(1, 20):
        print(i, ":", countDashes(a, i)+countDashes(b, i))

