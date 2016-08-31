#!/usr/bin/env bash
set -u #exit on using uninitialised vars
set -o pipefail # exit if error somewhere in pipe
set -e #exit on execution error

# paths
SCRIPTPATH=$1
CPUCORES=$3
REFNAME=$2
NODEPATH=$4

DATAPATH=$NODEPATH/Data
REFPATH=$NODEPATH/References
DESTPATH=$NODEPATH/Analysis
DORICPATH=$NODEPATH/DoriC

SUBJOB=5
SUBCPU=$[$CPUCORES/$SUBJOB]

echo "Running alignment stage"
cd $DATAPATH

# pick up .gz, .bz2 or .fastq files and prune file ending
declare -a files

shopt -s nullglob
fileEndings=(".fastq.bz2" ".fastq.gz" ".fastq")

for ending in "${fileEndings[@]}"; do
    files=( *_1"$ending" )
    files=( "${files[@]%_1$ending}" )
    fileEnding="$ending"
    if [ ! ${#files[@]} -eq 0 ]; then
        break
    fi
done

for fn in "${files[@]}"; do
	bowtie2 -k 15 --score-min L,-0.6,-0.6 -p $CPUCORES -x $REFPATH/Index/$REFNAME -1 "$fn"_1"$fileEnding" -2 "$fn"_2"$fileEnding" -S "$fn".sam
done

parallel -j $CPUCORES rm -f {}_1"$fileEnding" {}_2"$fileEnding" ::: "${files[@]}"
parallel "mkdir {} && mv {}.sam {}" ::: "${files[@]}"

echo "Moving into subjobs"
parallel -j $SUBCPU $SCRIPTPATH/bin/buildHelper.sh $1 $2 $SUBCPU $4 {} ::: "${files[@]}"
