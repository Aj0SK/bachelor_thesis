# reference file that we use to generate signal
refFilePath = "../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../data/kmer_model.hdf5"

repeatSignal = 10

fromRef, toRef = 100000, 100050
contig = "contig1"

import sys
import numpy as np
from pyfaidx import Fasta
from nadavca.dtw import KmerModel

sys.path.append("../helpers/hypothesis")
from signalHelper import stringToSignal

import matplotlib.pyplot as plt

################################################################################
# load reference sequence and create index for fast mapping

ref = Fasta(refFilePath)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

sequence = str(ref[contig][fromRef:toRef])

signal = stringToSignal(sequence, mod, repeatSignal=10)

plt.plot(signal)
plt.show()
