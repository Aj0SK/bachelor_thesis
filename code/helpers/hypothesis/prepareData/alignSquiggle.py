refFilePath = "../../../data/sapIngB1.fa"
readsFilePath = "../../../data/pos-basecalled"
outFilePath = "./alignedSquiggles.txt"

maxReads = 500

import os
import sys
import nadavca
from pyfaidx import Fasta

sys.path.append("../")
from signalHelper import getReadsInFolder

reads = getReadsInFolder(readsFilePath, minSize=0)

outFile = open(outFilePath, "w")

for read in reads:
    print(f"Ostava este {maxReads}")

    if maxReads == 0:
        break
    
    nadavca_align = nadavca.align_signal(refFilePath, [read])
    
    if nadavca_align == None:
        print(200*"x")
        continue
    
    if nadavca_align[0] == None:
        print(200*"y")
        continue
    
    if len(nadavca_align) != 1:
        continue

    maxReads -= 1
    readName = os.path.basename(read)
    nadavca_align = nadavca_align[0]
    
    fromSignal, toSignal = nadavca_align[0].signal_range
    fromRef, toRef = nadavca_align[0].reference_range
    ctg = nadavca_align[0].contig_name
    strand = -1 if nadavca_align[0].reverse_complement else 1
    #refSeq = "".join(nadavca_align[0].reference_part)
    
    outFile.write(f"{readName} {ctg} {strand} {fromRef} {toRef} {fromSignal} {toSignal}\n")
    
    #if strand == 1:
    #    r = str(Fasta(refFilePath)[ctg][fromRef:toRef])
    #else:
    #    r = str(-Fasta(refFilePath)[ctg][fromRef:toRef])

outFile.close()
