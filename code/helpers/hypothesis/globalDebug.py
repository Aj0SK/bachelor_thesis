################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

import argparse

parser = argparse.ArgumentParser(description='Options usable in this module.')
parser.add_argument("-p", "--sampleRead", help="location of positive sample read", type=str, default=None)
parser.add_argument("-n", "--sampleFakeRead", help="location of positive sample read", type=str, default = None)
args = parser.parse_args()

# here we want to chose some part of signal that is not close to beg or end of read
refFilePath = "../../data/sapIngB1.fa"
sampleRead = args.sampleRead
sampleFakeRead = args.sampleFakeRead

kmerModelFilePath = "../../data/kmer_model.hdf5"
# number of levels we use
numLevels = 6
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
from signalHelper import getLevels, stringToSignal, getSignalFromRead, getSeqfromRead, stringAllignment

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

def globalNormalize(signal):
    newSignal = copy.deepcopy(signal)
    for i in range(len(signal)):
        winBeg = max(0, i-winSize//2)
        winEnd = min(len(signal), i+winSize//2)
        newSignal[i] -= np.median(signal[winBeg:winEnd])
        newSignal[i] /= np.std(signal[winBeg:winEnd]-np.median(signal[winBeg:winEnd]), dtype="float64") + 0.000000001
        newSignal[i] = min(newSignal[i], maxSignal)
        newSignal[i] = max(newSignal[i], minSignal)
    return newSignal

############s####################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

################################################################################
# load basecalled sequence and signal
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastqSeq, basecallTable = getSeqfromRead(sampleRead)

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

hits = [i for i in sequenceIndex.map(fastqSeq) if i.q_en - i.q_st > 40*len(fastqSeq)//100]
print(hits[0].q_st)
print(hits[0].q_en)
print(len(fastqSeq))
if len(hits) != 1:
    print(len(hits))
    print("Bad number of hits")
    exit(0)
signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])
print("Of total len {0}, matched is [{1} {2}]".format(len(fastqSeq), hits[0].q_st, hits[0].q_en))

artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)

# load fake signal from read
fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = np.array(fakeSignal, dtype = float)

artifSignal = artifSignal[:15000]
originalSignal = originalSignal[:15000]
fakeSignal = fakeSignal[:15000]

artifSignal = globalNormalize(artifSignal)
originalSignal = globalNormalize(originalSignal)
fakeSignal = globalNormalize(fakeSignal)

artifString = "".join(getLevels(artifSignal,
                                kernelLen = kernelLen,
                                winSize = winSize,
                                numLevels = numLevels,
                                minSignal = minSignal,
                                maxSignal = maxSignal,
                                shift = winSize, normalize = False))
originalString = "".join(getLevels(originalSignal,
                                   kernelLen = kernelLen,
                                   winSize = winSize,
                                   numLevels = numLevels,
                                   minSignal = minSignal,
                                   maxSignal = maxSignal,
                                   shift = winSize, normalize = False))
fakeString = "".join(getLevels(fakeSignal,
                               kernelLen = kernelLen,
                               winSize = winSize,
                               numLevels = numLevels,
                               minSignal = minSignal,
                               maxSignal = maxSignal,
                               shift = winSize, normalize = False))

allig1, allig2, = stringAllignment(originalString, artifString)
allig3, allig4, = stringAllignment(fakeString, artifString)

print("\n{0} vs {1}".format(len(allig1), len(allig3)))

print(allig1[:500])
print(allig2[:500])
