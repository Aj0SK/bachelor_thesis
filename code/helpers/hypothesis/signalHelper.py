# Default settings, can be changed in the future as they can be passed to functions
# using them

# len of levelStrings
defaultKernelLen = 10
# len of vertical cut
defaultWinSize = 20 * defaultKernelLen
# number of horizontal cuts
defaultNumberOfLevels = 12
# minimal and maximal normalized signal we want to process
minSignal, maxSignal = -2.0, 2.0
# how long does one nucleotid passes through pore
repeatSignal = 8
# how many times we want to plot graph of normalized signal
showGraph = 0

import copy
import numpy as np
from itertools import groupby
import nadavca
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import h5py

# return basecalling info from read
def getSeqfromRead(filename):
    sequenceFile = h5py.File(filename, 'r')
    seq = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]
    basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]
    return seq, basecallEventTable

# return signal as list from single read
def getSignalFromRead(filename):
    readFile = h5py.File(filename, 'r')
    readName = str(list(readFile['Raw/Reads'])[0])
    rawData = readFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    return rawData

# Convert string in *ref_str* into signal(list of floats) using kmer_model loaded in *mod*
def stringToSignal(ref_str, mod):
    num_ref = [nadavca.alphabet.inv_alphabet[base] for base in ref_str]
    signal = mod.get_expected_signal(num_ref, [0]*4, [0]*4)
    signal = [x for x in signal for i in range(repeatSignal)]
    return signal

def normalizeWindow(w):
    w -= np.median(w)
    w /= np.std(w, dtype="float64") + 0.0000000001
    w[w<minSignal] = minSignal
    w[w>maxSignal] = maxSignal
    return

# cuts signal vertically into windows, normalizes windows and then cuts them
# horizontally into levels
def getLevels(signal, kernelLen = defaultKernelLen, winSize = defaultWinSize, numLevels = defaultNumberOfLevels):
    results = []
    
    # cut into windows and for every windows do the normalization and horiz. cutting
    for winBeg in range(0, signal.shape[0]-winSize, 1):
        winEnd = winBeg + winSize

        # normalize window
        currWindow = copy.deepcopy(signal[winBeg:winEnd])
        normalizeWindow(currWindow)

        levelSize = (maxSignal-minSignal)/numLevels + 0.0000000001

        # cut into horizontal levels
        outString = ""
        for i in currWindow:
            myLevel = int((i-minSignal)/levelSize)
            outString += chr(ord('a')+myLevel)

        outString = "".join([k for k, g in groupby(outString)])
        if len(outString) < kernelLen:
            print("* ", end = " ")
            continue
        # in case we want to graph our normalized signal
        '''
        global showGraph
        if showGraph is not 0:
            print(outString)
            plt.plot(currWindow)
            plt.show()
            showGraph -= 1
        '''
        results.append(outString[:kernelLen])
    return results

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
