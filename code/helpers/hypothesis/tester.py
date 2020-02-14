sampleRead = "../../data/albacore-output-pos/magnu_20181010_FAH93149_MN26672_sequencing_run_sapIng_19842_read_1000_ch_43_strand.fast5"
kmerModelFilePath = "../../data/kmer_model.hdf5"

import h5py
from signalHelper import getLevels, stringToSignal, normalizeWindow
from nadavca.dtw import KmerModel
import numpy as np
import matplotlib.pyplot as plt

def getRead(filename):
    readFile = h5py.File(filename, 'r')
    readName = str(list(readFile['Raw/Reads'])[0])
    rawData = readFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    return rawData

def getSignal(filename):
    sequenceFile = h5py.File(filename, 'r')
    seq = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]
    return seq

################################################################################
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

# read fastq from file
fastq = getSignal(sampleRead)
# read second line and decode it into string from bytes
fastqSeq = fastq.splitlines()[1].decode()
# create artificial signal
fabricatedSignal = np.array(stringToSignal(fastqSeq, mod), float)

# load original signal from read
originalSignal = getRead(sampleRead)
originalSignal = np.array(originalSignal, dtype = float)

# signal to list of level strings
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

for w in range(10, 20):
    fabrSignalLev = getLevels(fabricatedSignal, kernelLen = 10, winSize = w * 10)
    counter = 0
    d = {}
    for i in fabrSignalLev:
        if i in d:
            counter += 1
    print("\n{0} {1}".format(w, counter))
