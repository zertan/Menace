#!/usr/bin/env bash
#shopt -s nullglob
#files=( RefSeq/Fasta/*.fasta )
#files=( "${files[@]%_1.fastq.bz2}" )


# add ti numbers to references
while IFS=$'\t' read -r -a arr
do
	sed -i "" "s/^>gi/>ti|${arr[1]}|gi/" "${REF_PATH}/Fasta/${arr[0]}.fasta"
done < ${REF_PATH}/taxIDs.txt

#for fn in "${files[@]}"; do 
#	sed -i '' 's/>gi\|>ti\|/'$asd'\gi\|/g'
