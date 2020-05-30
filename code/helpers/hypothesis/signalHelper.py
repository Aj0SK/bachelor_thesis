# Default settings, can be changed in the future as they can be passed to functions
# using them

# len of levelStrings
defaultKernelLen = 7
# len of vertical cut
defaultWinSize = 15 * defaultKernelLen
# number of horizontal cuts
defaultNumberOfLevels = 5
# minimal and maximal normalized signal we want to process
defaultMinSignal, defaultMaxSignal = -2.0, 2.0
# how long does one nucleotid passes through pore
defaultRepeatSignal = 8
# how many times we want to plot graph of normalized signal
showGraph = 0

import copy
import numpy as np
from itertools import groupby
import nadavca
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import h5py

import os
import glob

# get reads from folder with size at least minSize(to ensure longer reads)
def getReadsInFolder(path, minSize=1000000):
    fileNames = glob.glob(path + "/**/*.fast5", recursive=True)
    fileNames = [i for i in fileNames if os.path.getsize(i) > minSize]
    return fileNames

def getAlignedIndex(path):
    index = {}
    with open(path, "r") as f:
        for line in f:
            line = [l.strip() for l in line.split()]
            readName, ctg, strand, fromRef, toRef, fromSignal, toSignal = line
            index[readName] = ctg, int(strand), int(fromRef), int(toRef), int(fromSignal), int(toSignal) 
    return index

# return basecalling info from read
def getSeqfromRead(filename):
    sequenceFile = h5py.File(filename, "r")
    seq = sequenceFile["/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"][()]
    basecallEventTable = sequenceFile[
        "/Analyses/Basecall_1D_000/BaseCalled_template/Events"
    ][()]
    seq = seq.decode("utf-8").split("\n")[1]
    return seq, basecallEventTable


# return nucleotide sequence that corresponds to signal from fromSignal to toSignal
def seqSignalCor(fromSignal, toSignal, basecallTable):
    signalFrTo = ""
    for i in Table_Iterator(basecallTable):
        if i[1] >= fromSignal and i[1] <= toSignal:
            signalFrTo += str(chr(i[0]))
    return signalFrTo


# return signal as list from single read
def getSignalFromRead(filename):
    readFile = h5py.File(filename, "r")
    readName = str(list(readFile["Raw/Reads"])[0])
    rawData = readFile["Raw/Reads/" + readName + "/" + "Signal"][()]
    return rawData


# Convert string in *ref_str* into signal(list of floats) using kmer_model loaded in *mod*
def stringToSignal(ref_str, mod, repeatSignal=defaultRepeatSignal):
    num_ref = [nadavca.alphabet.inv_alphabet[base] for base in ref_str]
    signal = mod.get_expected_signal(num_ref, [0] * 4, [0] * 4)
    signal = [x for x in signal for i in range(repeatSignal)]
    return signal


def normalizeWindow(
    w, minSignal=defaultMinSignal, maxSignal=defaultMaxSignal, shift=0, noise=False
):
    w -= np.median(w)
    w /= np.std(w, dtype="float64") + 0.0000000001
    # w /= np.median(np.abs(w)) + 0.0000000001
    # if noise == True: w += np.random.normal(0, 0.20 , w.shape[0])
    w += shift
    w[w < minSignal] = minSignal
    w[w > maxSignal] = maxSignal
    # w = medfilt(w, kernel_size = 7)
    # w = kf.smooth(w)[0].flatten()
    return

'''
def getLevelString(
    w,
    minSignal=defaultMinSignal,
    maxSignal=defaultMaxSignal,
    numLevels=defaultNumberOfLevels,
):
    outString = ""
    levelSize = (maxSignal - minSignal) / numLevels + 0.0000000001
    for i in w:
        myLevel = int((i - minSignal) / levelSize)
        outString += chr(ord("a") + myLevel)

    outString = "".join([k for k, g in groupby(outString)])
    return outString
'''

"""
def getLevelString(w, minSignal = defaultMinSignal, maxSignal = defaultMaxSignal, numLevels = defaultNumberOfLevels, borderArea=defaultBorderArea):
    outString = ""
    levelSize = (maxSignal-minSignal)/numLevels
    borders = [i for i in np.arange(minSignal+levelSize, maxSignal, levelSize)]
    for i in w:
        for j in range(len(borders)):
            if abs(borders[j]-i) < borderArea:
                outString += chr(ord('a')+numLevels+j)
                break
        else:
            myLevel = int((i-minSignal)/levelSize)
            outString += chr(ord('a')+myLevel)
            continue

    outString = "".join([k for k, g in groupby(outString)])
    return outString
"""

# cuts signal vertically into windows, normalizes windows and then cuts them
# horizontally into levels
def getLevels(
    signal,
    kernelLen=defaultKernelLen,
    winSize=defaultWinSize,
    numLevels=defaultNumberOfLevels,
    minSignal=defaultMinSignal,
    maxSignal=defaultMaxSignal,
    shift=1,
    normalize=True,
):
    results = []

    # cut into windows and for every windows do the normalization and horiz. cutting
    for winBeg in range(0, signal.shape[0] - winSize, shift):
        winEnd = winBeg + winSize

        # normalize window
        if normalize:
            currWindow = copy.deepcopy(signal[winBeg:winEnd])
            normalizeWindow(currWindow, minSignal, maxSignal)
        else:
            currWindow = signal[winBeg:winEnd]

        # cut into horizontal levels
        outString = getLevelString(currWindow, minSignal, maxSignal, numLevels)

        # in case we want to graph our normalized signal
        """
        global showGraph
        if showGraph is not 0:
            print(outString)
            plt.plot(currWindow)
            plt.show()
            showGraph -= 1
        """
        if len(outString) < kernelLen:
            # print("*", end = " ")
            continue
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
                self.localindex = 5 - int(self.table[self.tableindex][5])
            else:
                raise StopIteration
        self.totalindex += 1
        self.localindex += 1
        return (
            self.table[self.tableindex][4][self.localindex - 1],
            self.table[self.tableindex][1],
            self.table[self.tableindex][1],
        )


################################################################################
# allign two strings in O(n*m)
def stringAllignment(str1, str2):
    n, m = len(str1) + 1, len(str2) + 1
    tab = [[0] * m for _ in range(n)]

    for i in range(1, n):
        for j in range(1, m):
            if str1[i - 1] == str2[j - 1]:
                tab[i][j] = tab[i - 1][j - 1] + 1
            tab[i][j] = max(tab[i][j], tab[i][j - 1])
            tab[i][j] = max(tab[i][j], tab[i - 1][j])

    out1, out2 = "", ""
    i, j = n - 1, m - 1
    while i != 0 or j != 0:
        if i > 0 and tab[i][j] == tab[i - 1][j]:
            out1 += str1[i - 1]
            out2 += "-"
            i -= 1
        elif j > 0 and tab[i][j] == tab[i][j - 1]:
            out1 += "-"
            out2 += str2[j - 1]
            j -= 1
        else:
            out1 += str1[i - 1]
            out2 += str2[j - 1]
            i -= 1
            j -= 1
    return (out1[::-1], out2[::-1])


def countDashes(alligned, consec):
    alligned = "Z" + alligned + "Z"
    counter = 0
    expected = "-" * consec
    for i in range(0, len(alligned) - (consec + 2) + 1, 1):
        if (
            alligned[i] == "-"
            or alligned[i + consec + 1] == "-"
            or alligned[i + 1 : i + 1 + consec] != expected
        ):
            continue
        counter += 1
    alligned = alligned[1:-1]
    return counter


################################################################################

def getLevelString(signal, smoothParam, levels, overflow):
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

def computeNorm(signal, start, end):
    med = np.median(signal[start:end])
    std = np.std(signal[start:end] - med)
    return med, std


def computeString(signal, start, end, shift, scale, levels, overflow=0):
    min = -2 * scale + shift
    max = 2 * scale + shift
    step = (max - min) / levels

    lev = [min]
    for i in range(1, levels + 1):
        lev.append(lev[i - 1] + step)

    lev[0] = -10000
    lev[levels] = 10000

    outString = ""
    lastlevel = -1
    for s in signal:
        if (
            (lastlevel < 0)
            or (s < lev[lastlevel - 1] - overflow * step)
            or (s > lev[lastlevel] + overflow * step)
        ):
            level = 1
            while s > lev[level]:
                level += 1
            outString += chr(ord("a") + level - 1)
            lastlevel = level
    return outString


def computeStringInformative(
    signal, start, end, shift, scale, levels, overflow=0, repeatSignal=10
):
    min = -2 * scale + shift
    max = 2 * scale + shift
    step = (max - min) / levels

    lev = [min]
    for i in range(1, levels + 1):
        lev.append(lev[i - 1] + step)

    lev[0] = -10000
    lev[levels] = 10000

    outString = ""
    lastlevel = -1
    info = []
    for i in range(len(signal)):
        s = signal[i]
        if (
            (lastlevel < 0)
            or (s < lev[lastlevel - 1] - overflow * step)
            or (s > lev[lastlevel] + overflow * step)
        ):
            level = 1
            while s > lev[level]:
                level += 1
            outString += chr(ord("a") + level - 1)
            info.append(i // repeatSignal)
            lastlevel = level

    return outString, info


def smoothSignal(signal, window_len):
    newsignal = []
    sum = window_len * signal[0]
    for i in range(0, len(signal)):
        if i < window_len:
            sum -= signal[0]
        else:
            sum -= signal[i - window_len]
        sum += signal[i]
        newsignal.append(sum / window_len)

    return newsignal

def smoothSignalMed(signal, window_len):
    newsignal = []
    win = [signal[0]] * window_len
    for i in range(0, len(signal)):
        win = win[1:]
        win.append(signal[i])
        newsignal.append(np.median(win))

    return newsignal

def buildDictionary(string, k):
    dict = {}
    for i in range(0, len(string) - k + 1):
        kmer = string[i : i + k]
        if kmer not in dict:
            dict[kmer] = 0
        dict[kmer] += 1

    return dict


def overlappingKmers(string1, string2, k):
    dict1 = buildDictionary(string1, k)
    dict2 = buildDictionary(string2, k)
    intersect = 0
    for kmer in dict1:
        if kmer in dict2:
            intersect += 1
    return intersect


##################################

from random import randint


def produceRandom(size, levels):
    last = -1
    outStr = ""
    for i in range(size):
        g = last
        while g == last:
            g = randint(0, levels - 1)
        outStr += chr(ord("a") + g)
        last = g
    return outStr
