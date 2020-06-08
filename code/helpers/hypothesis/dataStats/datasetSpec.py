readsFilePath = "../../../data/pos-basecalled"

maxReads = 1000

import os
import sys
import numpy as np

sys.path.append("../")
from signalHelper import getReadsInFolder, getSignalFromRead, getSeqfromRead

reads = getReadsInFolder(readsFilePath, minSize=0)

signalLengths, seqLengths = [], []

for read in reads:
    if maxReads == 0:
        break

    try:
        signal = getSignalFromRead(read)
        seq, _ = getSeqfromRead(read)
    except:
        continue

    maxReads -= 1
    signalLengths.append(len(signal))
    seqLengths.append(len(seq))


meanSignal = np.mean(signalLengths)
medianSignal = np.median(signalLengths)

meanSeq = np.mean(seqLengths)
medianSeq = np.median(seqLengths)

print(f"Not found {maxReads} reads!")
print(f"Mean signal len is {meanSignal}")
print(f"Median signal len is {medianSignal}")
print(f"Mean sequence len is {meanSeq}")
print(f"Median sequence len is {medianSeq}")
