#!/usr/bin/env bash
#shopt -s nullglob
#files=( RefSeq/Fasta/*.fasta )
#files=( "${files[@]%_1.fastq.bz2}" )


# add ti numbers to references
while IFS=$'\t' read -r -a arr
do
	LINE="$(head -1 "$1/Fasta/${arr[0]}.fasta")"
	
	if [[ $(uname) == 'Linux' ]]; then
		LINE2=$(echo "$LINE" | sed "s/^>/>ti|${arr[1]}|org|/")
		echo "$LINE2"
		if [ "$LINE"!="$LINE2" ]; then
			( echo "$LINE2" && tail -n +2 "$1/Fasta/${arr[0]}.fasta" ) > file.new && mv -f file.new "$1/Fasta/${arr[0]}.fasta"
		fi
			#sed -i "1s/.*/${LINE}/" "$1/Fasta/${arr[0]}.fasta"
	elif [[ $(uname) == 'Darwin' ]]; then
		sed -i "" "s/^>/>ti|${arr[1]}|org|/" "$1/Fasta/${arr[0]}.fasta"
	fi
	#fi
done < $1/taxIDs.txt
