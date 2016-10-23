#!/usr/bin/env bash
set -u #exit on using uninitialised vars
set -o pipefail # exit if error somewhere in pipe
set -e #exit on execution error

function mvAll {
	mv $1 $(echo $1 | sed -E 's/.+(N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+).+/\1\.bam/')
}
export -f mvAll

# paths
SCRIPTPATH=$1
CPUCORES=$3
REFNAME=$2
NODEPATH=$4
FILE=$5
SUBJOBS=$6

DATAPATH=$NODEPATH/Data
REFPATH=$NODEPATH/References
DESTPATH=$NODEPATH/Analysis
DORICPATH=$NODEPATH/DoriC

function getcpu {
	pathojobs=$(pgrep "pathoscope ID" | wc -l)
	if [[ $SUBJOBS == $pathojobs ]]; then
		echo $CPUCORES
	else
		endcores=$[$pathojobs*($CPUCORES-1)/($SUBJOBS-$pathojobs)]
		endcores=$((2*$CPUCORES<$endcores?2*$CPUCORES:$endcores))
		endcores=$(($endcores<2?$CPUCORES:$endcores))
		echo $endcores
	fi
}

cd $FILE

echo "$FILE: Running pathoscope stage"
pathoscope ID -alignFile $FILE.sam -expTag $FILE -fileType sam
rm -f $FILE.sam

echo "$FILE: Running coverage stage"
samtools view -@ $CPUCORES -bS updated_$FILE.sam -o $FILE.bam
rm -f updated_$FILE.sam

echo "$FILE: Sorting"
samtools sort -@ $CPUCORES $FILE.bam -o $FILE.sorted.bam
rm -f $FILE.bam

bamtools split -refPrefix $FILE -in $FILE.sorted.bam -reference
rm -f $FILE.sorted.bam

#parallel -j $(getcpu) "mv {.}* '{= s:\.[^.]+$::;s:\.[^.]+$::; =}'" ::: *.bam

parallel mvAll {} ::: *.bam

echo "$FILE: Calculating coverage"
perl -v > /dev/null
rc=$?;
if [[ $rc != 0 ]]; then
	parallel samtools mpileup -q8 {} '>' {.}.d ::: *.bam
	parallel "cat {} | awk '{print \$2,\$4}' > {}epth"  ::: *.d
	rm -f *.d
else
	parallel -j 1 samtools view -@ $CPUCORES -q8 -o {.}.sam {} ::: *.bam
	parallel $SCRIPTPATH/bin/interp.pl {} $REFPATH/Headers/{.}.xml '>' {.}.depth ::: *.sam
	rm -f *.sam
fi
rm -f *.bam

find . -size -100 -type f -name \*.depth

echo "$FILE: Performing piecewise fits"

parallel $SCRIPTPATH/bin/piecewiseFit.py {} $REFPATH/Headers/ $DORICPATH/ ::: *.depth

endings=("log" "png" "npy" "depth")
for ending in ${endings[@]}; do
	mkdir $ending && files=(*.$ending) && if [ "${#files[@]}" -ne 0 ]; then mv *.$ending $ending; fi
done

cd depth

echo "$FILE: Performing secondary fits"
cat $REFPATH/taxIDs.txt | awk 'BEGIN {FS=" "} {print $2}' | sort -u | parallel "\$SCRIPTPATH/bin/test.py \$REFPATH {} \$(grep -r \" {}\$\" \$REFPATH/taxIDs2.txt |  awk 'BEGIN {FS=\" \"} {print \$1}' | tr '\n' ' ')"

parallel "\$SCRIPTPATH/bin/piecewiseFit.py {} \$REFPATH/Headers/ \$DORICPATH/" ::: *.npy

cd ..

mkdir log2 && mv depth/*.log log2
mkdir png2 && mv depth/*.png png2
mkdir npy2 && mv depth/*.npy npy2
