#!/usr/bin/env bash
set -u #exit on using uninitialised vars
set -o pipefail # exit if error somewhere in pipe
set -e #exit on execution error

# paths
REFNAME=$2
NODEPATH=$4
DATAPATH=$NODEPATH/Data
REFPATH=$NODEPATH/$REFNAME
DESTPATH=$NODEPATH/Analysis
DORICPATH=$NODEPATH/DoriC
SCRIPTPATH=$1
CPUCORES=$3

echo "Running alignment stage"
cd $DATAPATH

# pick up .gz or .bz2 files and prune file ending
declare -a files

shopt -s nullglob
fileEndings=(".fastq.bz2" ".fastq.gz" ".fastq")

for ending in "${fileEndings[@]}"; do
    files=( *"$ending" )
    files=( "${files[@]%_1$ending}" )
    fileEnding="$ending"
    if [ ! ${#files[@]} -eq 0 ]; then
        break
    fi
done

for fn in "${files[@]}"; do
	bowtie2 -k 5 -p $CPUCORES -x $REFPATH/Index/$REFNAME -1 "$fn"_1"$fileEnding" -2 "$fn"_2"$fileEnding" -S "$fn".sam
done

parallel rm -f {}_1.fastq.bz2 {}_2.fastq.bz2 ::: "${files[@]}"

echo "Running pathoscope"

parallel mkdir {} ::: "${files[@]}"

parallel pathoscope ID -alignFile {}.sam -expTag {} -fileType sam -outDir {} ::: "${files[@]}"
rm -f *.sam

for RUNNAME in "${files[@]}"; do
	cd $RUNNAME
	echo "$RUNNAME: Running coverage stage"

    samtools view -@ $CPUCORES -bS updated_$RUNNAME.sam -o $RUNNAME.bam
    rm -f updated_$RUNNAME.sam

    echo "$RUNNAME: Sorting"
    samtools sort -@ $CPUCORES $RUNNAME.bam $RUNNAME.sorted
    rm -f $RUNNAME.bam

    bamtools split -in $RUNNAME.sorted.bam -reference
    mv $RUNNAME.sorted.bam ../

    echo "$RUNNAME: Calculating coverage"
    parallel samtools mpileup -q14 {} '>' {.}.depth ::: *.bam
    rm -f *.bam

    echo "$RUNNAME: Performing prefit work"

    parallel "cat {} | awk '{print \$2,\$4}' > {}2"  ::: *.depth

    rm -f *.depth
    find . -size 0 -delete

    for f in $(ls --color=none | grep .depth2); do
        mv $f $(echo $f |  sed -r 's/.+(N[CT]_[0-9]{6}\.[0-9]+).+/\1\.depth/')
    done

    echo "$RUNNAME: Performing piecewise fits"

    parallel $SCRIPTPATH/bin/piecewiseFit.py {} $REFPATH/Headers/ $DORICPATH/ ::: *.depth

    mkdir log
    mv *.log log
	
	mkdir png
	mv *.png png
	
	mkdir npy
	mv *.npy npy

	mkdir depth
	mv *.depth depth

	cd $DATAPATH
done

