kernelLen = 10
winSize = 15 * kernelLen
refFilePath = "../../data/sapIngB1.fa"
kmerModelFilePath = "../../data/kmer_model.hdf5"
minSignal, maxSignal = -2.0, 2.0
defaultNumberOfLevels = 12
repeatSignal = 8
showGraph = 0

################################################################################

import sys
import numpy as np
from pyfaidx import Fasta
import nadavca
from nadavca.dtw import KmerModel
from itertools import groupby
import copy
import matplotlib.pyplot as plt

def stringToSignal(ref_str, mod):
    num_ref = [nadavca.alphabet.inv_alphabet[base] for base in ref_str]
    signal = mod.get_expected_signal(num_ref, [0]*4, [0]*4)
    signal = [x for x in signal for i in range(repeatSignal)]
    return signal

def getLevels(signal, numLevels = defaultNumberOfLevels):
    results = []

    for winBeg in range(0, signal.shape[0]-winSize, winSize):
        winEnd = winBeg + winSize

        currWindow = copy.deepcopy(signal[winBeg:winEnd])
        currWindow -= np.median(currWindow)
        currWindow /= np.std(currWindow, dtype="float64") + 0.0000000001
        currWindow[currWindow<minSignal] = minSignal
        currWindow[currWindow>maxSignal] = maxSignal

        #minSignal, maxSignal = np.amin(currWindow), np.amax(currWindow)
        levelSize = (maxSignal-minSignal)/numLevels + 0.0000000001

        outString = ""

        for i in currWindow:
            myLevel = int((i-minSignal)/levelSize)
            outString += chr(ord('a')+myLevel)

        outString = "".join([k for k, g in groupby(outString)])
        if len(outString) < kernelLen:
            print("*", end = " ")
            continue

        global showGraph
        if showGraph is not 0:
            print(outString)
            plt.plot(currWindow)
            plt.show()
            showGraph -= 1

        results.append(outString[:kernelLen])
    return results

################################################################################
genes = Fasta(refFilePath)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

genesSignal = {}

for i in genes:
    genesSignal[i.name] = stringToSignal(str(i), mod)

hashTable = {}

for key in genesSignal.keys():
    out = getLevels(np.array(genesSignal[key], dtype = float))
    for i in out:
        if i not in hashTable:
            hashTable[i] = 0
        hashTable[i] += 1

print("Hashtable ready!")
########################################
### hashTable

import h5py
import glob

def getRead(filename):
    readFile = h5py.File(filename, 'r')
    readName = str(list(readFile['Raw/Reads'])[0])
    rawData = readFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    return rawData

def countMatch(hTable, signalStrings):
    counter = 0
    for i in signalStrings:
        if i in hTable:
            counter += 1
    return counter

########################################

readsPosFilePath = "../../data/pos"
readsNegFilePath = "../../data/neg"

posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

print("Positive:")

for filePath in posFast5:
    posRead = np.array(getRead(filePath), dtype = float)
    oneReadLevels = getLevels(posRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(negRead)), end=' ')

print()
print("Negative:")

for filePath in negFast5:
    negRead = np.array(getRead(filePath), dtype = float)
    oneReadLevels = getLevels(negRead)
    print(countMatch(hashTable, oneReadLevels), end='')
    print("/" + str(len(negRead)), end=' ')

print()
