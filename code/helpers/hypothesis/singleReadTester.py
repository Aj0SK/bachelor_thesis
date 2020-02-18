sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"
kmerModelFilePath = "../../data/kmer_model.hdf5"

import h5py
from signalHelper import getLevels, stringToSignal, normalizeWindow, getSignalFromRead, getSeqfromRead
from nadavca.dtw import KmerModel
import numpy as np
import math

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

class Table_Iterator:
    def __init__(self, basecallEventTable):
        self.table = basecallEventTable
        self.tableindex = 0
        self.localindex = 0
        self.totalindex = 0

    def __iter__(self):
        return self

    def __next__(self):
        while self.localindex == 5:
            if self.tableindex + 1 != len(self.table):
                self.tableindex += 1
                self.localindex = 5-int(self.table[self.tableindex][5])
            else:
                raise StopIteration
        self.totalindex += 1
        self.localindex += 1
        return self.table[self.tableindex][4][self.localindex-1], self.table[self.tableindex][1], self.table[self.tableindex][1]

################################################################################
# load signal and basecalled sequence
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastq, basecallTable = getSeqfromRead(sampleRead)
# read second line and decode it into string from bytes
fastqSeq = fastq.splitlines()[1].decode()
# create artificial signal
fabricatedSignal = np.array(stringToSignal(fastqSeq, mod), float)

# load original signal from read
originalSignal = getSignalFromRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

################################################################################
# get only part of signal and find corresponding basecalled sequence
x = Table_Iterator(basecallTable)

fr, to = 24000, 24500
signalFrTo = ""

for i in x:
    #print(chr(i[0]) + " " + str(i[1]) + " " + str(i[2]))
    if i[1] >= fr and i[1] <= to:
        signalFrTo += str(chr(i[0]))

fabricatedSignal = np.array(stringToSignal(signalFrTo, mod), float)
originalSignal = originalSignal[fr:to]

normalizeWindow(originalSignal)
normalizeWindow(fabricatedSignal)

################################################################################
# modify signal
from scipy.signal import medfilt

betterSignal = np.array(getSignalFromRead(sampleRead), dtype = float)[fr:to]
normalizeWindow(betterSignal)

betterSignal = medfilt(betterSignal, kernel_size = 9)

################################################################################
# graph original, fabricated and modified signal
o = np.linspace(0, originalSignal.shape[0], originalSignal.shape[0])
f = np.linspace(0, fabricatedSignal.shape[0], fabricatedSignal.shape[0])
e = np.linspace(0, betterSignal.shape[0], betterSignal.shape[0])

fig, axs = plt.subplots(3)
fig.suptitle('Original, edited original and fabricated')
axs[0].plot(o, originalSignal, label='data 1', color = 'r')
axs[1].plot(e, betterSignal, label='data 2', color = 'g')
axs[2].plot(f, fabricatedSignal, label='data 3', color = 'b')

axs[0].hlines(y=[-2+i*(4/12) for i in range(12)], xmin=0, xmax=len(o), linewidth=1, color='gray')
axs[1].hlines(y=[-2+i*(4/12) for i in range(12)], xmin=0, xmax=len(e), linewidth=1, color='gray')
axs[2].hlines(y=[-2+i*(4/12) for i in range(12)], xmin=0, xmax=len(f), linewidth=1, color='gray')
plt.show()

levelO = getLevels(originalSignal)
levelF = getLevels(fabricatedSignal)
levelE = getLevels(betterSignal)

print(levelO)
print(levelF)
print(levelE)

print(signalFrTo)
# signal to list of level strings
'''
origSignalLev = getLevels(originalSignal)
fabrSignalLev = getLevels(fabricatedSignal)

d = {}

for i in origSignalLev:
    d[i] = 1

counter = 0
for i in fabrSignalLev:
    if i in d:
        counter += 1

print("{0}/{1}".format(counter, len(fabrSignalLev)))
'''
