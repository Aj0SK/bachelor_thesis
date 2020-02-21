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
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

ref = Fasta(refFilePath)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

################################################################################

import sys
import numpy as np
import glob
from signalHelper import getSeqfromRead, getSignalFromRead, getLevels, Table_Iterator, stringToSignal

# load filenames of all positive and negative reads
posFast5 = glob.glob(readsPosFilePath + '/*.fast5', recursive=True)
negFast5 = glob.glob(readsNegFilePath + '/*.fast5', recursive=True)

for filePath in posFast5:
    seq, basecallTable = getSeqfromRead(filePath)
    x = Table_Iterator(basecallTable)
    originalSignal = np.array(getSignalFromRead(filePath), dtype = float)
    oneReadLevel = getLevels(originalSignal)
    
    hits = sequenceIndex.map(seq.splitlines()[1].decode())
    for hit in hits:
        refString = str(ref[hit.ctg][hit.r_st:hit.r_en])
        signal = stringToSignal(refString, mod)
        hitLevels = getLevels(np.array(signal, dtype = float))
    
    print(len([1 for i in hitLevels if i in oneReadLevel]))
    print("Done " + filePath)
