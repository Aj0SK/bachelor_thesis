minimalHitLength = 30
seqLimit = 1000
refFilePath = "../../data/sapIngB1.fa"
readsFilePath = "../../data/reads.fastq"
kmerModelFilePath = "../../data/kmer_model.hdf5"

import sys
import numpy as np
import mappy as mp
from pyfaidx import Fasta
import h5py
import nadavca
from nadavca.dtw import KmerModel

sequenceIndex = mp.Aligner(refFilePath)
assert sequenceIndex, "failed to load/build reference index"

hits = []

with open(readsFilePath, "r") as f: 
    lines = f.readlines()
    old = [lines[i+1] for i in range(len(lines)) if lines[i][0] == '@']
    counter = 0
    for seq in old:
        if counter >= seqLimit:
            break
        refHits = [i.ctg for i in list(sequenceIndex.map(seq))]
        if len(set(refHits)) is not 1:
            continue
        hits.append((refHits[0], seq))
        counter += 1

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

for contig, ref_str in hits:
    if len(ref_str) < minimalHitLength:
        continue
    num_ref = [nadavca.alphabet.inv_alphabet[base] for base in ref_str if base is not '\n']
    signal = mod.get_expected_signal(num_ref, [0, 0, 0, 0], [0, 0, 0, 0])
    print(contig)
    print(signal)



