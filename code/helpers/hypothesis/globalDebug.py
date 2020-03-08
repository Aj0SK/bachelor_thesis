################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

# here we want to chose some part of signal that is not close to beg or end of read
fromSignal, toSignal = 19000, 21000
fromSignalFake, toSignalFake = 20000, 22000
refFilePath = "../../data/sapIngB1.fa"
sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"
sampleFakeRead = "../../data/albacore-output-neg/magnu_20181218_FAH93149_MN26672_sequencing_run_sapFun_DNase_flush_88184_read_1001_ch_42_strand.fast5"

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
from signalHelper import getGlobalLevels, stringToSignal, getSignalFromRead, getSeqfromRead, stringAllignment

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.close('all')

def globalNormalize(signal):
    newSignal = copy.deepcopy(signal)
    for i in range(len(signal)):
        winBeg = max(0, i-winSize)
        winEnd = min(len(signal), i+winSize)
        newSignal[i] -= np.median(signal[winBeg:winEnd])
        newSignal[i] /= np.std(signal[winBeg:winEnd]-np.median(signal[winBeg:winEnd]), dtype="float64")
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
fastq, basecallTable = getSeqfromRead(sampleRead)
# read second line and decode it into string from bytes
fastqSeq = fastq.splitlines()[1].decode()

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

hits = list(sequenceIndex.map(fastqSeq))
assert len(hits) == 1, "Too many hits"
signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])
print("Of total len {0} matched is [{1} {2}]".format(len(fastqSeq), hits[0].q_st, hits[0].q_en))

artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)

# load original signal from read
fakeSignal = getSignalFromRead(sampleFakeRead)
fakeSignal = np.array(fakeSignal, dtype = float)

#artifSignal = artifSignal[:15000]
#originalSignal = originalSignal[:15000]
#fakeSignal = fakeSignal[:15000]

artifSignal = globalNormalize(artifSignal)
originalSignal = globalNormalize(originalSignal)
fakeSignal = globalNormalize(fakeSignal)

artifString = "".join(getGlobalLevels(artifSignal, kernelLen = kernelLen, winSize = winSize, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal, shift=winSize))
originalString = "".join(getGlobalLevels(originalSignal, kernelLen = kernelLen, winSize = winSize, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal, shift=winSize))
fakeString = "".join(getGlobalLevels(fakeSignal, kernelLen = kernelLen, winSize = winSize, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal, shift=winSize))

allig1, allig2, = stringAllignment(originalString, artifString)
allig3, allig4, = stringAllignment(fakeString, artifString)

print("\n{0} vs {1}".format(len(allig1), len(allig3)))
