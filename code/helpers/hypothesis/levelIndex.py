################################################################################
# This module takes one basecalled read, then choses some part of signal and finds
# corresponding basecalled sequence. This sequence is then alligned to reference
# sequence to get rid of basecalling errors. Now when we have real signal and
# corresponding reference sequence we can try to create fake signal and see the
# differences.

refFilePath = "../../data/sapIngB1.fa"
kmerModelFilePath = "../../data/kmer_model.hdf5"

repeatSignal = 10
overflow = 0.3

levels = range(4, 15)

import numpy as np

from pyfaidx import Fasta
from nadavca.dtw import KmerModel
from signalHelper import stringToSignal
from signalHelper import smoothSignal, computeNorm, computeString

ref = Fasta(refFilePath)
mod = KmerModel.load_from_hdf5(kmerModelFilePath)

for contig in ref:
    refSeqPos = str(contig[:])
    refSeqNeg = str(contig[:].complement)

    refSignalPos = np.array(
        stringToSignal(refSeqPos, mod, repeatSignal=repeatSignal), float
    )
    refSignalNeg = np.array(
        stringToSignal(refSeqNeg, mod, repeatSignal=repeatSignal), float
    )

    refSignalPos = smoothSignal(refSignalPos, 5)
    refSignalNeg = smoothSignal(refSignalNeg, 5)

    refSignalPosShift, refSignalPosScale = computeNorm(
        refSignalPos, 0, len(refSignalPos)
    )
    refSignalNegShift, refSignalNegScale = computeNorm(
        refSignalNeg, 0, len(refSignalNeg)
    )

    for l in levels:
        refStringPos = computeString(
            refSignalPos,
            0,
            len(refSignalPos),
            refSignalPosShift,
            refSignalPosScale,
            l,
            overflow=overflow,
        )

        refStringNeg = computeString(
            refSignalNeg,
            0,
            len(refSignalNeg),
            refSignalNegShift,
            refSignalNegScale,
            l,
            overflow=overflow,
        )

        print(f"{contig.name} {l} + {refStringPos}")
        print(f"{contig.name} {l} - {refStringNeg}")
