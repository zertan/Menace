#!/usr/bin/env bash
set -u #exit on using uninitialised vars
set -o pipefail # exit if error somewhere in pipe
set -e #exit on execution error

# paths
REFNAME=gemQin
DATAPATH=$TMPDIR/Data
REFPATH=$TMPDIR/RefSeq/$REFNAME
DESTPATH=$TMPDIR/Analysis
DORICPATH=$TMPDIR/DoriC
SCRIPTPATH=$SLURM_SUBMIT_DIR

CPUCORES=16
THREADMEMORY=1800

echo "Running alignment stage"
cd $DATAPATH

shopt -s nullglob
files=( *_1.fastq.bz2 )
files=( "${files[@]%_1.fastq.bz2}" )

parallel bzip2 -d {} ::: *.bz2

for fn in "${files[@]}"; do
	gem-mapper -p -b -I $REFPATH/Index/$REFNAME.gem -T $CPUCORES -q offset-33 --gem-quality-threshold 26 -s 0 -d all -D 1 -1 "$fn"_1.fastq -2 "$fn"_2.fastq -o "$fn"
	sleep 10
	gem-2-sam --max-memory 1000000 -l -I $REFPATH/Index/$REFNAME.gem -q offset-33 -i "$fn".map -o "$fn".sam
	#gemtools convert -I $REFPATH/Index/$REFNAME.gem -i "$fn".map -o "$fn".bam -q 33 -t $CPUCORES --no-index --no-sort -p
	rm -f "$fn".map
	#samtools view -@ $CPUCORES "$fn".bam -o "$fn".sam
	#rm -f "$fn".bam
done

parallel rm -f {}_1.fastq {}_2.fastq ::: "${files[@]}"

parallel mkdir {} ::: "${files[@]}"

echo "Running pathoscope"

parallel pathoscope ID -alignFile {}.sam -expTag {} -fileType sam ::: "${files[@]}"
parallel rm -f {}.sam ::: "${files[@]}"

for RUNNAME in "${files[@]}"; do
	echo "$RUNNAME: Running coverage stage"

    samtools view -@ $CPUCORES -F 4 -bS updated_$RUNNAME.sam -o $RUNNAME.bam
    rm -f updated_$RUNNAME.sam

    echo "$RUNNAME: Sorting"
    samtools sort -@ $CPUCORES $RUNNAME.bam $RUNNAME.sorted
    rm -f $RUNNAME.bam
done

parallel "samtools view -H {}.sorted.bam > {}.header" :::  "${files[@]}"

parallel "cat {}.header | sed -r 's/\///g' > {}.reheader" :::  "${files[@]}"

parallel "samtools reheader {}.reheader {}.sorted.bam > {}.sorted2.bam" ::: "${files[@]}"

parallel bamtools split -in {}.sorted2.bam -reference ::: "${files[@]}"
rm -f *.sorted2.bam *.sorted.bam
rm -f *.header *.reheader

parallel "rsync -avm --remove-source-files --include='{}*' -f 'hide,! */' ./     ./{}" ::: "${files[@]}"

awkfun() {
awk 'BEGIN { FS = "\t" } ; {print $2,$4}'
}

for RUNNAME in "${files[@]}"; do
	cd $RUNNAME

    echo "$RUNNAME: Calculating coverage"
    parallel samtools mpileup -q10 {} '>' {.}.depth ::: *.bam
    rm -f *.bam

    echo "$RUNNAME: Performing prefit work"

    parallel "cat {} | awk -F $'\t' '{print \$2,\$4}' > {}2"  ::: *.depth

    rm -f *.depth
    find . -size 0 -delete
	
 	files2=( *.depth2 )

    for f in "${files2[@]}"; do
        mv "$f" $(echo "$f" |  sed -r 's/.+(N[CT]_[0-9]{6}\.[0-9]+).+/\1\.depth/')
    done
	
	parallel "cat {} | awk '{sum+=\$2} END { print \"{.}\",sum}'" ::: *.depth > coverage.csv
	
    echo "$RUNNAME: Performing piecewise fits"

    $SCRIPTPATH/refSel.py ./ $REFPATH/headers/

	parallel $SCRIPTPATH/piecewiseFit.py {}.depth $REFPATH/headers/ $DORICPATH/ :::: referenceACC.txt

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

