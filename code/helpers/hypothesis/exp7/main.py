# index global level string of reference and try to find read level string in it
# using minimap

# reference file that we use to generate signal
refFilePath = "../../../data/sapIngB1.fa"
refIndex = "../../../data/ref.index"
kmerModelFilePath = "../../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

posTestCases, negTestCases = 100, 0
level = 4
repeatSignal = 10
overflow = 0.3
smoothParam = 5
fromRead, toRead = 5000, 10000
contigNum = 5

levelStringFa = "helper392478942.fa"

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

from signalHelper import stringAllignment, countDashes

def getGlobalString(signal):
    signal = smoothSignal(signal, smoothParam)
    signalShift, signalScale = computeNorm(signal, 0, len(signal))
    signalString = computeString(signal,
                                 0,
                                 len(signal),
                                 signalShift,
                                 signalScale,
                                 level,
                                 overflow=overflow)
    return signalString


################################################################################

refFasta = Fasta(refFilePath)

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
lengths = {}

if os.path.exists(levelStringFa):
    os.remove(levelStringFa)

with open(refIndex, "r") as outFile:
    for line in outFile:
        line = line.split()
        if int(line[1]) != level or line[2] != "+":
            continue
        levelStr = line[3]
        contigName = line[0]
        #storeContig[line[0]] = levelStr
        processed.append(contigName)
        helperDict = {"a": "A", "b": "C", "c": "G", "d": "T"}
        refString = "".join(helperDict[i] for i in levelStr)
        lengths[contigName] = len(refString)
        with open(levelStringFa, "a") as f:
            f.write(">" + contigName + "\n")
            f.write(refString.upper() + "\n")
        contigNum -= 1
        if contigNum == 0:
            break

index = mp.Aligner(levelStringFa, w = 25, best_n = 5, extra_flags = 0x100000)
#index = mp.Aligner(levelStringFa, w = 25)
assert index, "failed to load/build reference index"

print("Hashtable readyfor {0} nums!".format(contigNum))
print(f"Processed: {str(processed)}")
#######################################

def processRead(path, contig_name = None, refPosition = -1):
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
        print(f"Strand is {hit.strand}")
        diff = (hit.r_en-hit.r_st)-(hit.q_en-hit.q_st)
        print("{0}: {1} vs {2}".format(hit.ctg, hit.r_en-hit.r_st, hit.q_en-hit.q_st))
        hitPosition = hit.r_st/lengths[hit.ctg]
        print(f"Position of hit is {hitPosition}")
        if abs(hitPosition-refPosition) <= 0.01:
            if contig_name != None and hit.ctg != contig_name:
                print("Zle urceny contig!")
                print("Return False.")
                #return False
            return True
        '''if diff < 0.05*(hit.q_en-hit.q_st):
            a, b = stringAllignment(str(refFasta[hit.ctg][hit.r_st:hit.r_en]), helperString[hit.q_st:hit.q_en])
            for i in range(1, 20):
                print(i, ":", countDashes(a, i)+countDashes(b, i))

            if contig_name != None and hit.ctg != contig_name:
                print("Zle urceny contig!")
                print("Return False.")
                return False
            print("Return True.")
            return True'''
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
    
    refPosition = hits[0].r_st/len(refFasta[hits[0].ctg])
    print(f"Position in reference is {refPosition}")
    if processRead(filePath, contig_name = hits[0].ctg, refPosition = refPosition) == True:
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
