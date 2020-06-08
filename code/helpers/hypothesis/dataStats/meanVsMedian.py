refFilePath = "../../../data/sapIngB1.fa"

readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

readCount = 100

import sys
import mappy as mp
import numpy as np
from pyfaidx import Fasta

sys.path.append("../")
from signalHelper import getReadsInFolder, getSignalFromRead


import matplotlib
import matplotlib.pyplot as plt

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

posReads = getReadsInFolder(readsPosFilePath, minSize=0)
negReads = getReadsInFolder(readsNegFilePath, minSize=0)

################################################################################

totalCount = 0

suma = 0

for readFile in posReads:
    if totalCount == readCount:
        break

    try:
        readSignal = getSignalFromRead(readFile)
    except:
        continue

    mean = np.mean(readSignal)
    median = np.median(readSignal)
    stdev = np.std(readSignal)

    totalCount += 1

    # print(f"Mean is {mean}")
    # print(f"Median is {median}")
    #print(f"Stdeviation is {stdev}")

    suma += abs(mean - median)

    #print(f"Diff is {mean-median}")

print(f"Average diff is {suma/totalCount}")
