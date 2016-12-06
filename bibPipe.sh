#!/usr/bin/env bash

# on docker img

#### fix headers
# first strip referenes from : and - in header then:

### check minimal LCB length -c

########## run in docker img
#mugsy -p mult -c 300 Fasta/*.fasta #NC_000913.3.fasta NC_007779.1.fasta
#cp /tmp/mult.maf .
###########

# extract seqs of equal length from maf and fix headers
# maf file handling: bx-python


#maf_thread_for_species.py
#maf_limit_to_species.py $(maf_species_in_all_files.py mult.maf) < mult.maf > mult_lim.maf

#maf_to_fasta.py < mult_lim.maf > mult.fasta

rm -f mult_lim.fasta
maf_to_concat_fasta.py $(maf_species_in_all_files.py mult.maf) < mult.maf > mult_lim.fasta

#strip first set, do for each set
#head -6 mult_lim.fasta > contained.fasta

# run Gblocks for each set
Gblocks mult_lim.fasta -t=d
# build tree from concatenated seqs

#FastTree -nt contained.fasta-gb > tree
rm -f tree
FastTree -nt -gtr -nosupport < mult_lim.fasta-gb > tree

# choose which refs to drop
rm -f tree_pruned
refs=$(~/Documents/GitRepos/PTRloc/PTR-Pipeline/bin/prune_ref.py tree) #| tr '\n' ','
# refs=$(~/Documents/GitRepos/PTRloc/PTR-Pipeline/bin/prune_ref.py tree2)

## run program to get refs under certain phylgen distance shown to be distinguishable

## drop refs from multiple alignment file  mult.fasta

# 


#######
# from generalized refs, create 'artifical map to coverage genome'


###### old with mauve ####################

#maf_to_fasta.py < /tmp/mult.maf >> other.fsa

############## 

#progressiveMauve --output=full_alignment.xmfa genome1.fasta genome2.fasta genome3.fasta genome4.fasta

#stripSubsetLCBs full_alignment.xmfa full_alignment.xmfa.bbcols core_alignment.xmfa 500 1

##########################################

# choose to concat or not
#maf_to_concat_fasta.py $(echo ${refs[*]} | tr ' ' ',') < (mult_lim.maf | sed '/^>/!s/-//g') > mult.fasta

#maf_limit_to_species.py $(echo ${refs[*]} | tr ' ' ',') < mult.maf | maf_to_fasta.py | sed '/^>/!s/-//g' > mult.fasta

## choose lcb limited to all or not 


#maf_to_fasta.py < mult.maf > mult.fasta

#sed '/^>/!s/-//g' mult.fasta > core_alignment_gapless.fasta
#sed 's/=//g' core_alignment_gapless_x.fasta > core_alignment_gapless_xx.fasta

#sed -r -i 's/^>.*:([0-9]+).*/>\1/g' core_alignment_gapless_xx.fasta

#bowtie2-build `parallel echo {}.1.fasta, ::: "${refs[@]}"` cag2

### make artificial sample
##

rm -f *.fq
rm -f out*
rm -f tot_*

parallel wgsim -1 100 -2 100 -N 1000 Fasta/{}.1.fasta {}_1.fq {}_2.fq ::: "${refs[@]}"

cat *_1.fq > tot_1.fq
cat *_2.fq > tot_2.fq

rm -f mult_cat.fasta

#cat `parallel echo {}.1.fasta ::: "${refs[@]}"` > mult_cat.fasta

cat `parallel echo {} ::: Fasta/*` > mult_cat.fasta

#wgsim -N 100000 NZ_HG738867.fasta out1.fq out2.fq
#wgsim -N 100000 NC_012892.fasta out3.fq out4.fq
#wgsim -N 100000 NC_011601.fasta out5.fq out6.fq

#cat out1.fq out3.fq out5.fq > tot_1.fq
#cat out2.fq out4.fq out6.fq > tot_2.fq

bowtie2 -x cag2 -1 tot_1.fq -2 tot_2.fq -a --very-sensitive-local -p 4 | samtools view -b - > output.bam

# do this for sub problems

# bamtools split -reference -in output.bam

# grep ">" core_alignment_gapless.fasta | sed -E 's/>(.*)/output.REF_\1.bam/' > split_names

# split -l 3 split_names 

# parallel bamtools merge -list {} -out {}.bam ::: x* # check with -p 1

# split -a 4 -p 'label=' mult_lim.maf 

#cat xaaaa xaaab | maf_to_fasta.py - | sed '/^>/!s/-//g' > xaa.fasta

#parallel -j 1 cat Fasta/{}.1.fasta ::: "${refs[@]}" > mult_from_refs.fasta

rm -f genome_info.tr
parseAlignment output.bam -o alignment_info.prob --trSeqFile mult_cat.fasta --trInfoFile genome_info.tr --uniform --verbose #-l [depend on transcripts] --mateNamesDiffer
estimateVBExpression -o abundance alignment_info.prob -t genome_info.tr --saveAlignmentProbs

# send to gurobi assignment

#pathoscope ID -alignFile alignment_file.sam -expTag alignment_file -fileType sam

#tail -n +6 updated_tot.sam | awk '{print $1 $3}' | sed -E 's/.*(NC_[0-9]{6}\.[0-9]+).*(\1).*/\1 \2/;tx;d;:x' | wc -l

#tail -n +6 updated_tot.sam | awk '{print $1 $3}' | grep -E '.*NC_002695\.1.*NC_002695\.1' | wc -l

# check inds

../gurobi_assign.py $(pwd)

#cat genome_info.tr-NEW | sed -E 's/none (.+)\.[A-Z].+/\1/g' 

# tail -n +2 genome_info.tr-NEW | sed -E 's/none (.+)\.[A-Z].+/\1/g' > col1

#./compareSim.pl abundance.m_phi | grep -e '.*NC_002695\.1.* 2' | wc -l


# short gets more!!!!!!!!

# docker run interactive with host mounted folder:
# docker run -it -v /Users/Shared:/mnt/share centos/hedani /bin/bash

# save docker state
# sudo docker ps -l
# Commit changes to the container:
# sudo docker commit <container_id> <container name> 


# FastTree
# 