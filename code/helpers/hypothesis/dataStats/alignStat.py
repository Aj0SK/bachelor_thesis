refFilePath = "../../../data/sapIngB1.fa"

readsPosFilePath = "../../../data/pos-basecalled"
readsNegFilePath = "../../../data/neg-basecalled"

ratios = [0.001, 0.10, 0.25, 0.5, 0.75, 0.95]
readCount = 500

import sys
import mappy as mp
import numpy as np
from pyfaidx import Fasta

sys.path.append("../")
from signalHelper import getSeqfromRead, getReadsInFolder


import matplotlib
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 18})

################################################################################

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

posReads = getReadsInFolder(readsPosFilePath, minSize=0)
negReads = getReadsInFolder(readsNegFilePath, minSize=0)

################################################################################

negHitsByRatios = [0] * len(ratios)
totalCount = 0

for readFile in negReads:
    if totalCount == readCount:
        break

    try:
        readFastq, _ = getSeqfromRead(readFile)
    except:
        continue

    hits = [aln for aln in referenceIdx.map(readFastq)]
    hits = [(hit.q_en - hit.q_st) for hit in hits]

    totalCount += 1

    for i in reversed(range(len(ratios))):
        hit = max(hits + [0])
        if hit >= ratios[i] * len(readFastq):
            negHitsByRatios[i] += 1
            # break

posHitsByRatios = [0] * len(ratios)
totalCount = 0

for readFile in posReads:
    if totalCount == readCount:
        break

    try:
        readFastq, _ = getSeqfromRead(readFile)
    except:
        continue

    hits = [aln for aln in referenceIdx.map(readFastq)]
    hits = [(hit.q_en - hit.q_st) for hit in hits]

    totalCount += 1

    for i in reversed(range(len(ratios))):
        hit = max(hits + [0])
        if hit >= ratios[i] * len(readFastq):
            posHitsByRatios[i] += 1
            # break

print(posHitsByRatios)
print(negHitsByRatios)

# assert sum(posHitsByRatios) == readCount
# assert sum(negHitsByRatios) == readCount

# taken from https://matplotlib.org/3.2.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py

labels = [">=" + str(i * 100) + "%" for i in ratios]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width / 2, posHitsByRatios, width, label="Saprochaete ingens reads")
rects2 = ax.bar(
    x + width / 2, negHitsByRatios, width, label="Saprochaete fungicola reads"
)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("Number of alignments")
ax.set_xlabel("Ratio of alignment length to read length")
ax.set_title("Alignment of reads to Saprochaete ingens reference")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate(
            "{}".format(height),
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(0, 3),  # 3 points vertical offset
            textcoords="offset points",
            ha="center",
            va="bottom",
        )


autolabel(rects1)
autolabel(rects2)

# fig.tight_layout()

plt.show()
