#!/usr/bin/env bash
fe="xml"
files=( *.$fe )
declare -A a

for f in "${files[@]}"; do
	gi=$(head -2 $f | grep '<Id>' | sed -r "s/.*<Id>(.+)<\/Id>.*/\1/")
	acc=$(cat $f | grep '<Item Name="Extra" Type="String">' | sed -r "s/.+(N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+).+/\1/")
	a[$gi]=$acc
	#mv $f "$acc.$fe"
	#mv $f $(head -1 $f | sed -r "s/.+(N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+).+/\1\.$fe/")
done

#for i in "${!a[@]}"
#do
#	echo "key  : $i"
#	if ["$i" -eq "814569384"]; then echo "match"; fi
#	echo "value: ${a[$i]}"
#done

cols=( $( awk '{print $1}' $1 ) )
clen=${#cols[@]}

declare -a final

for col in "${cols[@]}"; do
	#if test "${a[$col]+isset}"; then echo "yes"; else echo "no"; fi
	#if ["$col" -eq "$col"] 2>/dev/null; then
	#	echo number
	#fi
	#printf "%s" $col
	final+=(${a[$col]})
	#echo ${a[$col]}
	#printf "%s" "${a[$col]}"
done

#for (( i=0; i<${clen}; i++ ));
#do
#	tmp=${cols[$i]}
#	final[$i]=${a[$tmp]}
#	echo "--"
#	echo ${a[$tmp]}
#done

printf "%s\n" "${final[@]}"	> $1"2"

paste $1 $1"2" > $1"3"

#variable=$(printf "%s\n" "${final[@]}")

awk -F $'\t' 'BEGIN {OFS = FS} {print $4"\t"$2"\t"$3}' $1"3" > $1"4"

#sed -i 'Ns/.*/replacement-line/' file.txt


