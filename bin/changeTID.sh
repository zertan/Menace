#!/usr/bin/env bash
#shopt -s nullglob
#files=( RefSeq/Fasta/*.fasta )
#files=( "${files[@]%_1.fastq.bz2}" )


# add ti numbers to references
while IFS=$'\t' read -r -a arr
do
	#LINE=$(head -1 "$1/Fasta/${arr[0]}.fasta")
	#if [[ "$LINE" =~ "^>gi.+" ]]
	#then
	if [[ $(uname) == 'Linux' ]]; then
		sed -i "s/^>gi/>ti|${arr[1]}|gi/" "$1/Fasta/${arr[0]}.fasta"
	elif [[ $(uname) == 'Darwin' ]]; then
		sed -i "" "s/^>gi/>ti|${arr[1]}|gi/" "$1/Fasta/${arr[0]}.fasta"
	fi
	#fi
done < $1/taxIDs.txt

#for fn in "${files[@]}"; do 
#	sed -i '' 's/>gi\|>ti\|/'$asd'\gi\|/g'