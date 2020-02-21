################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

# here we want to chose some part of signal that is not close to beg or end of read
fromSignal, toSignal = 24000, 27000
refFilePath = "../../data/sapIngB1.fa"
sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"

kmerModelFilePath = "../../data/kmer_model.hdf5"
# number of levels we use
numLevels = 12
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 10

import math
import h5py
import numpy as np
import mappy as mp
from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import getLevels, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead, Table_Iterator

import matplotlib
# this allows us to see graphs through ssh
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

################################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

################################################################################
# load signal and basecalled sequence
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastq, basecallTable = getSeqfromRead(sampleRead)
# read second line and decode it into string from bytes
fastqSeq = fastq.splitlines()[1].decode()

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

################################################################################
# get only part of signal and find corresponding basecalled sequence, then find
# it in reference
signalFrTo = ""

for i in Table_Iterator(basecallTable):
    #print(chr(i[0]) + " " + str(i[1]) + " " + str(i[2]))
    if i[1] >= fromSignal and i[1] <= toSignal:
        signalFrTo += str(chr(i[0]))

# try to find string in reference
hits = list(sequenceIndex.map(signalFrTo))
assert len(hits) == 1, "Too many hits"
signalFrTo = str(ref[hits[0].ctg][hits[0].r_st:hits[0].r_en])

# create artificial signal
artifSignal = np.array(stringToSignal(signalFrTo, mod, repeatSignal = repeatSignal), float)
originalSignal = originalSignal[fromSignal:toSignal]

################################################################################
# modify original signal
from scipy.signal import medfilt

procSignal = np.array(getSignalFromRead(sampleRead), dtype = float)[fromSignal:toSignal]
procSignal = medfilt(procSignal, kernel_size = 9)

levelO = getLevels(originalSignal, numLevels = numLevels, minSignal = minSignal, maxSignal = maxSignal)
levelO = set(levelO)
levelE = getLevels(procSignal, minSignal = minSignal, maxSignal = maxSignal)
levelE = set(levelE)

levelF = getLevels(artifSignal, minSignal = minSignal, maxSignal = maxSignal)

# count have many strings are shared by fabricatedSignal and original/processedSignal
counterO = len([i for i in levelF if i in levelO])
counterE = len([i for i in levelF if i in levelE])

print("Artificial match rate with original signal {0}/{1}".format(counterO, len(levelF)))
print("Artificial match rate with preprocessed signal {0}/{1}".format(counterE, len(levelF)))
#print(levelO)
#print(levelE)
#print(levelF)
#print(signalFrTo)

################################################################################
# Graph normalized original, processed and artificial signal

normalizeWindow(originalSignal)
normalizeWindow(artifSignal)
normalizeWindow(procSignal)

o = np.linspace(0, originalSignal.shape[0], originalSignal.shape[0])
f = np.linspace(0, artifSignal.shape[0], artifSignal.shape[0])
e = np.linspace(0, procSignal.shape[0], procSignal.shape[0])

fig, axs = plt.subplots(3)
fig.suptitle('Original, modified original(processed) and artificial signal')
axs[0].plot(o, originalSignal, label='data 1', color = 'r')
axs[1].plot(e, procSignal, label='data 2', color = 'g')
axs[2].plot(f, artifSignal, label='data 3', color = 'b')

axs[0].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(o), linewidth=1, color='gray')
axs[1].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(e), linewidth=1, color='gray')
axs[2].hlines(y=[minSignal+i*((maxSignal-minSignal)/numLevels) for i in range(numLevels+1)], xmin=0, xmax=len(f), linewidth=1, color='gray')
plt.show()
