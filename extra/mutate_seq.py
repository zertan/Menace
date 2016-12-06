#!/usr/bin/env python
import numpy as np
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

def mutateSeq(mutable_seq,rate):
    nr=['A','C','G','T']
    nt={'A': 0,'C': 1,'G': 2,'T': 3}
    
    prob=np.matrix('0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0')/float(3)
    prob=prob.tolist()
    rand_ch=np.random.choice(range(len(mutable_seq)),np.round(len(mutable_seq)*rate))
    for i,c in enumerate(rand_ch):
        mutable_seq[c]=np.random.choice(nr,None,False,prob[nt[mutable_seq[c]]])
    
    return mutable_seq


handle = open(sys.argv[1], "rU")
for rec in SeqIO.parse(handle, "fasta"):
    seq = rec.seq
mutable_seq = MutableSeq(str(rec.seq), IUPAC.unambiguous_dna)
handle.close()

rec.seq=mutateSeq(mutable_seq,.01)

rec.id="mut_"+rec.id
rec.description="mut_"+rec.description
rec.name="mut_"+rec.name

with open(sys.argv[2], "w") as output_handle:
    SeqIO.write(rec, output_handle, "fasta")