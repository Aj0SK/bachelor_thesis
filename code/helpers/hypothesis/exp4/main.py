# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

posTestCases, negTestCases = 50, 50
levels = 4
repeatSignal = 10
kmerLength = 21
overflow = 0.3
smoothParam = 5
fromRead, toRead = 7000, 60000
contigNum = 2

################################################################################

import os
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

processed = []

if os.path.exists("helper.fa"):
    os.remove("helper.fa")

for contig in Fasta(refFilePath):
    print("Processing " + contig.name)
    processed.append(contig.name)
    contigSignal = stringToSignal(str(contig), mod, repeatSignal=repeatSignal)
    refString = getGlobalString(contigSignal)
    helperDict = {"a": "A", "b": "C", "c": "G", "d": "T"}
    refString = "".join(helperDict[i] for i in refString)
    with open("helper.fa", "a") as f:
        f.write(">" + contig.name + "\n")
        f.write(refString.upper() + "\n")
    contigNum -= 1
    if contigNum == 0:
        break

index = mp.Aligner("helper.fa")
assert index, "failed to load/build reference index"

print("Hashtable readyfor {0} nums!".format(contigNum))
#######################################

def processRead(path, contig_name = None):
    print(path)
    readSignal = np.array(getSignalFromRead(path), dtype=float)
    readSignal = readSignal[fromRead:toRead]
    readLevelString = getGlobalString(readSignal)
    
    helperDict = {"a": "A", "b": "C", "c": "G", "d": "T"}
    helperString = "".join(helperDict[i] for i in readLevelString)
    levelHits = list(index.map(helperString))
    
    if len(levelHits) == 0:
        print("Return False.")
        return False
    
    for hit in levelHits:
        diff = (hit.r_en-hit.r_st)-(hit.q_en-hit.q_st)
        print("{0}: {1} vs {2}".format(hit.ctg, hit.r_en-hit.r_st, hit.q_en-hit.q_st))
        if diff < 0.05*(hit.q_en-hit.q_st):
            print("Return True.")
            if contig_name != None and hit.ctg != contig_name:
                print("Zle urceny contig!")
            return True
    print("Return False.")
    return False

########################################

good, bad = 0, 0

print("Positive:")

counter = 0

for filePath in posFast5:
    #print(filePath)
    if counter == posTestCases:
        break
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
    
    counter += 1
    
    if processRead(filePath, contig_name = hits[0].ctg) == True:
        good += 1
    else:
        bad += 1
    print("{0}/{1}".format(good, good + bad))

print("\n\nNegative:")

counter = 0

for filePath in negFast5:
    #print(filePath)
    if counter == negTestCases:
        break
    try:
        readSeq, basecallTable = getSeqfromRead(filePath)
    except:
        continue
    if len(readSeq) < (toRead // repeatSignal):
        continue
    counter += 1
    if processRead(filePath) == False:
        good += 1
    else:
        bad += 1
    print("{0}/{1}".format(good, good + bad))

print("{0}/{1}".format(good, good + bad))
