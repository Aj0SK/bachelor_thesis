# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../data/kmer_model.hdf5"

# positive and negative reads folder
readsPosFilePath = "../../data/albacore-output-pos"
readsNegFilePath = "../../data/albacore-output-neg"

################################################################################
# build index on reference
import mappy as mp

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

################################################################################

import sys
import numpy as np
import glob
from signalHelper import getSeqfromRead

# load filenames of all positive and negative reads
posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

for filePath in posFast5:
    seq, table = getSeqfromRead(filePath)
    hits = sequenceIndex.map(seq.splitlines()[1].decode())
    for hit in hits:
        print("{} {}\t{}\t{}".format(hit.ctg, filePath, hit.r_st, hit.r_en))
