# Try to evaluate positive and negative reads based on the match count with our
# builded table

# reference file that we use to generate signal
refFilePath = "../../../data/sapIngB1.fa"
# kmer model
kmerModelFilePath = "../../../data/kmer_model.hdf5"

readsFolder = "/projects2/nanopore-2016/saping-20181010/reads/20181010_0936_sapIng/fast5/"

readsFastq = "/projects2/nanopore-2016/saping-20181010/reads.fastq"

maxReads = 50

targetContig = 'contig1'
targetBeg, targetEnd = 10000, 20000

import glob
import mappy as mp
from nadavca.dtw import KmerModel

################################################################################

def overlap(b, e, c, d):
    if e <= c or b >= d:
        return False
    return True

referenceIdx = mp.Aligner(refFilePath)
assert referenceIdx, "failed to load/build reference index"

mod = KmerModel.load_from_hdf5(kmerModelFilePath)

reads = glob.glob(readsFolder + "**/*.fast5", recursive=True)

with open(readsFastq, "r") as f:
    infoLine = False
    info = None
    counter = 0
    for line in f.readlines():
        if counter > maxReads:
            break
        if line[0] == "@":
            infoLine = True
            info = line
            #print(line)
            continue
        if infoLine == True:
            infoLine = False
            #print(line[:20])
            #print(line[-20:])
            hits = list(referenceIdx.map(line))
            goodHits = [i for i in hits if i.ctg == targetContig and i.strand == 1 and overlap(targetBeg, targetEnd, i.r_st, i.r_en)]
            if len(goodHits) != 0:
                read = [i[5:] for i in info.split() if i[:4] == "read"]
                ch = [i[3:] for i in info.split() if i[:2] == "ch"]
                #print(ch[0])
                #print(read)
                matchingFileNames = [i for i in reads if i.find("read_"+read[0] + "_ch_"+ch[0]) != -1]
                if len(matchingFileNames) == 1:
                    print(matchingFileNames[0])
                    counter += 1

            info = None
