kernelLen = 10
winSize = 15 * kernelLen
defaultNumberOfLevels = 12
minSignal, maxSignal = -2.0, 2.0
repeatSignal = 8
showGraph = 0

import copy
import numpy as np
from itertools import groupby
import nadavca
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

        levelSize = (maxSignal-minSignal)/numLevels + 0.0000000001

        outString = ""

        for i in currWindow:
            myLevel = int((i-minSignal)/levelSize)
            outString += chr(ord('a')+myLevel)

        outString = "".join([k for k, g in groupby(outString)])
        if len(outString) < kernelLen:
            print("*", end = " ")
            continue
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
