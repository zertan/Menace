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

if [ "$fileEnding" == ".fastq.bz2" ]; then
	parallel bzip2 -d {} ::: *.bz2
	rm -f *.bz2
fi

if [ "$fileEnding" == ".fastq.gz" ]; then
	parallel gzip -d {} ::: *.gz
	rm -f *.gz
fi

for fn in "${files[@]}"; do
	gem-mapper -p -b -I $REFPATH/Index/$REFNAME.gem -T $CPUCORES -q offset-33 --gem-quality-threshold 26 -s 0 -d all -D 1 -1 "$fn"_1.fastq -2 "$fn"_2.fastq -o "$fn"
	sleep 3
	gem-2-sam --max-memory 100 -l -I $REFPATH/Index/$REFNAME.gem -q offset-33 -i "$fn".map -o "$fn".sam
	#gem-2-sam --max-memory 100 -q offset-33 -i "$fn".map -o "$fn".sam
	#gemtools convert -I $REFPATH/Index/$REFNAME.gem -i "$fn".map -o "$fn".bam -q 33 -t $CPUCORES --no-index --no-sort -p
	rm -f "$fn".map
	#samtools view -@ $CPUCORES "$fn".bam -o "$fn".sam
	#rm -f "$fn".bam
done

parallel rm -f {}_1.fastq {}_2.fastq ::: "${files[@]}"

parallel mkdir {} ::: "${files[@]}"

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

	files=( *.bam )
	for f in "${files[@]}"; do
		mv $f $(echo $f | sed -r 's/.+(N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+).+/\1\.bam/')
	done

    echo "$RUNNAME: Calculating coverage"
	perl -v
    rc=$?; 
	if [[ $rc != 0 ]]; then 
		parallel samtools mpileup -q8 {} '>' {.}.d ::: *.bam
		parallel "cat {} | awk '{print \$2,\$4}' > {}epth"  ::: *.d
		rm -f *.d
	else
		parallel samtools view -@ $CPUCORES -q8 -o {.}.sam {} ::: *.bam
		parallel $SCRIPTPATH/bin/interp.pl {} $REFPATH/Headers/{.}.xml '>' {.}.depth ::: *.sam
		rm -f *.sam
	fi
	rm -f *.bam

    #echo "$RUNNAME: Performing prefit work"

    #parallel "cat {} | awk '{print \$2,\$4}' > {}2"  ::: *.depth

    #rm -f *.depth
    find . -size -100 -type f -name \*.depth 

    #for f in $(ls --color=none | grep .depth2); do
    #    mv $f $(echo $f |  sed -r 's/.+(N[CT]_[0-9]{6}\.[0-9]+).+/\1\.depth/')
    #done

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
