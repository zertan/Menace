#!/usr/bin/env bash

#IN1=comm00_1.fastq
#IN2=comm00_2.fastq

#IN1=sub1.fq
#IN2=sub2.fq

python PTRC.py DB -metadata meta -fasta multi.fasta -out comm0/sim_db

python PTRC.py -db_path_name comm0/sim_db CA -pe -i1 comm0/Data/$1 -i2 comm0/Data/$2 -m comm0/Data/sm.map -outfol comm0/

python PTRC.py -db_path_name comm0/sim_db PTR -infol comm0/ -o comm0/Out/ptrc_out -minpred 1 -csv_output -preds predetermined

#seqtk sample -s100 ERR525696_1.fastq 200000 > comm0/Data/sub1.fq
#seqtk sample -s100 ERR525696_2.fastq 200000 > comm0/Data/sub2.fq

#GTTGCGCGCCAATATCGATCAGGTCAAGCTGAACCGCAAGGTCAATGCGCTGGTGCGTGACGTGGATCTCGGTCTCGATATCGAGGATCTGACGTTCGGC
#CCCACTGCTCGAACTGTGCCGCGGATTCAATAACCGTGATTTTCGGTGTGTCCAAACCGGAACCCAGCTCGTCTTCTGCAGGACCATCAGACTGTGCGGC

#gem-mapper -I comm0/sim_db.gem -1 comm0/Data/sub1.fq -2 comm0/Data/sub2.fq -o comm0/Data/sm -q offset-33 -d all -p --max-extendable-matches all --max-extensions-per-match 5 -v

#gem-mapper -I comm0/sim_db.gem -1 comm0/Data/sub01.fq -2 comm0/Data/sub02.fq -o comm0/Data/sm -q offset-33 -d all -p --max-extendable-matches all --max-extensions-per-match 5 -v
