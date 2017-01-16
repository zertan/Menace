#!/usr/bin/env bash
set -u #exit on using uninitialised vars
set -o pipefail # exit if error somewhere in pipe
set -e #exit on execution error

# paths
SCRIPTPATH=$1
CPUCORES=$3
REFNAME=$2
NODEPATH=$4
OUTPATH=$5

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
fileEndings=(".fastq.gz" ".fastq.bz2" ".fastq" ".fq.bz2" ".fq.gz" ".fq")

for ending in "${fileEndings[@]}"; do
    files=( *_1"$ending" )
    files=( "${files[@]%_1$ending}" )
    fileEnding="$ending"
    if [ ! ${#files[@]} -eq 0 ]; then
        break
    fi
done

for fn in "${files[@]}"; do
	echo "Aligning $fn:"
	>&2 echo "$fn:"
	bowtie2 -k 4 --very-sensitive -p $CPUCORES -x $REFPATH/Index/$REFNAME -1 "$fn"_1"$fileEnding" -2 "$fn"_2"$fileEnding" -S "$fn".sam #--un-gz single_un_"$fn".fastq.gz --un-conc-gz paired_un_"$fn".fastq.gz --no-unal
done

#parallel -j $CPUCORES mv single_un_{}.fastq.gz $OUTPATH ::: "${files[@]}"
#parallel -j $CPUCORES mv paired_un_{}.fastq.1.gz $OUTPATH ::: "${files[@]}"
#parallel -j $CPUCORES mv paired_un_{}.fastq.2.gz $OUTPATH ::: "${files[@]}"

parallel -j $CPUCORES rm -f {}_1"$fileEnding" {}_2"$fileEnding" ::: "${files[@]}"
parallel "mkdir {} && mv {}.sam {}" ::: "${files[@]}"

echo "Moving into subjobs"
parallel -j $SUBJOB $SCRIPTPATH/bin/buildHelper.sh $1 $2 $SUBCPU $4 {} $SUBJOB ::: "${files[@]}"

echo "Finnished. Data in" $OUTPATH
