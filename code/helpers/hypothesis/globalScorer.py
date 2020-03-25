import sys
import numpy as np
import mappy as mp
from nadavca.dtw import KmerModel
from pyfaidx import Fasta
from signalHelper import getSeqfromRead, seqSignalCor, getSignalFromRead, stringToSignal, getReadsInFolder

def computeNorm(signal,start,end):
    med = np.mean(signal[start:end])
    std = np.std(signal[start:end]-med)
    return med,std

def computeString(signal,start,end,shift,scale,levels,overflow=0):
    min = -2*scale + shift
    max = 2*scale + shift
    step = (max-min) / levels

    lev = [min]
    for i in range(1,levels+1):
        lev.append(lev[i-1] + step)

    lev[0] = -10000
    lev[levels] = 10000
    
    outString = ""
    lastlevel = -1
    for s in signal:
        if (lastlevel < 0) or (s < lev[lastlevel-1]-overflow*step) or (s > lev[lastlevel]+overflow*step):
            level = 1
            while (s > lev[level]):
                level+=1
            outString += chr(ord('a')+level-1)
            lastlevel = level            
    return outString


def smoothSignal(signal,window_len):
    newsignal = []
    sum = window_len * signal[0]
    for i in range(0,len(signal)):
        if (i<window_len):
            sum -= signal[0]
        else:
            sum -= signal[i-window_len]
        sum += signal[i]
        newsignal.append(sum/window_len)

    return newsignal;

def buildDictionary(string,k):
    dict = {}
    for i in range(0,len(string)-k+1):
        kmer = string[i:i+k]
        if kmer not in dict:
            dict[kmer] = 0
        dict[kmer] += 1

    return dict;

def overlappingKmers(string1,string2,k):
    dict1 = buildDictionary(string1,k)
    dict2 = buildDictionary(string2,k)
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
    return intersect
            
    

##################################

refFile = "../../data/sapIngB1.fa"
readsPosFilePath = "../../data/pos-basecalled"
kmerModelFilePath = "../../data/kmer_model.hdf5"

maxTests = 150

levels = 5
repeatSignal = 8

signalFrom = 6000#int(sys.argv[3])
signalTo = 10000#int(sys.argv[4])

mod = KmerModel.load_from_hdf5(kmerModelFilePath)
posReads = getReadsInFolder(readsPosFilePath, minSize = 0)

totalG = [0]*9
totalF = [0]*9
badClas = [0]*9

print(len(posReads))

for readFile in posReads[:min(len(posReads), maxTests)]:
    print(readFile)
    ### read read
    readFastq, readEvents = getSeqfromRead(readFile)
    readSeq = seqSignalCor(signalFrom, signalTo, readEvents)
    readSignal = np.array(getSignalFromRead(readFile)[signalFrom:signalTo],
                             dtype=float)
    # readSeq - sequence cut out from read
    # readSignal - corresponding signal section

    ### get corresponding part of the reference using minimap2
    referenceIdx = mp.Aligner(refFile)
    assert referenceIdx, "failed to load/build reference index"
    # require a single hit with at least 95% coverage of length
    hits = [aln for aln in referenceIdx.map(readSeq)
            if aln.q_en - aln.q_st > 0.95*len(readSeq)]
    if len(hits) != 1:
        print("Too many or too few hits, skipping read.")
        continue
    hit = hits[0]

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
    
    print(readSeq)
    print(refSeq)
    print(fakeSeq)
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
    for i in [3,5,7,9,11,13,15,17,19]:
        x = overlappingKmers(readString2Sm,refString2Sm,i)
        totalG[((i-1)//2)-1] += x
        print(i,":",x,end="\t")
    print("")

    print("read vs. fake")
    for i in [3,5,7,9,11,13,15,17,19]:
        x = overlappingKmers(readString2Sm,fakeString2Sm,i)
        totalF[((i-1)//2)-1] += x
        print(i,":",x,end="\t")
    print("")
    
    print(totalG)
    print(totalF)
    for i in [3,5,7,9,11,13,15,17,19]:
        x = overlappingKmers(readString2Sm,refString2Sm,i)
        y = overlappingKmers(readString2Sm,fakeString2Sm,i)
        if x<y:
            badClas[((i-1)//2)-1] += 1
    
    print([totalG[i]/totalF[i] for i in range(len(totalG))])
    print(badClas)
