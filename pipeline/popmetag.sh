#!/bin/bash
#
# SNPs from metagenomic samples
#
# ./popmetag.sh [superfile] [reference] [name]
#
# - superfile: path of station/sample corresponding file
# - reference: fasta file of the reference
# - name: output name
#-------------------------------------------------------------------------------


# Declare and initialize
DIR=./mapping
TRIMMOMATIC=./bioinfo_tools/Trimmomatic-0.33/trimmomatic-0.33.jar
ANALYSIS=genpop_$3'_'$(date "+%Y-%m-%d")
RAW=$DIR/raw_data_fastq
ASSEMBLY=$2
REF=$DIR/reference_assembly/$ASSEMBLY

FILTER_ID=$DIR/filter95.py
FILTER_LCBAD=$DIR/filter_lc-bad.py
FILTER_LCGOOD=$DIR/filter_lc-good.py

# Create analysis repositories
mkdir $ANALYSIS
if [ ! -d filter_fastq ] ; then
  mkdir filter_fastq
fi
mkdir $ANALYSIS/alignment
mkdir $ANALYSIS/alignment_filt
mkdir $ANALYSIS/variants/
mkdir $ANALYSIS/variants/ref

# Create an index for assembly fasta file
cp $REF $ANALYSIS/variants/ref
samtools faidx $ANALYSIS/variants/ref/$ASSEMBLY

# Create coverage summary file
echo "STATION SAMPLE TOT H_SUM V_SUM H_COV V_COV" > $ANALYSIS/coverage_$3

# Number of stations in the analysis
N=$(cat $1 | wc -l)
echo "Total number of metagenomic samples to process: "$N

for i in $(seq 1 $N)
do
	STATION=$(sed -n $i'p' $1 | cut -f1)
	ECH=$(sed -n $i'p' $1 | cut -f2)

	echo "Processing station: "$STATION

## Quality Control  ------------
	# Trimmed bad quality extremities
	if [ ! -f $DIR/filter_fastq/$STATION'_'$SAMPLE* ] ; then
		echo "1) TRIMMING: "$STATION
		java -jar $TRIMMOMATIC PE -threads 30 -phred33 $DIR/raw_data_fastq/$SAMPLE'_1.fastq.gz'  $DIR/raw_data_fastq/$SAMPLE'_2.fastq.gz' $DIR/filter_fastq/$STATION'_'$SAMPLE'_1-trimmed.fastq.gz' $DIR/filter_fastq/$STATION'_'$SAMPLE'_output_forward_unpaired.fq.gz'  $DIR/filter_fastq/$STATION'_'$SAMPLE'_2-trimmed.fastq.gz' $DIR/filter_fastq/$STATION'_'output_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 >  $DIR/filter_fastq/'log_'$STATION'_'$SAMPLE'_trimmomatic.txt' 2>&1
	fi

# Alignment ------------
	echo "2) ALIGNMENT : "$STATION
	mkdir $ANALYSIS/alignment/$STATION

	bwa index $REF
	bwa mem -M -t 30 $REF $DIR/filter_fastq/$STATION'_'$SAMPLE'_1-trimmed.fastq.gz' $DIR/filter_fastq/$STATION'_'$SAMPLE'_2-trimmed.fastq.gz' > $ANALYSIS/alignment/$STATION/$SAMPLE.sam

	samtools view -Shu -F 4 $ANALYSIS/alignment/$STATION/$SAMPLE.sam | samtools sort - > $ANALYSIS/alignment/$STATION/$SAMPLE.sort.bam
	samtools index $ANALYSIS/alignment/$STATION/$SAMPLE.sort.bam
	rm $ANALYSIS/alignment/$STATION/$SAMPLE.sam

## Filtering ------------
	mkdir $ANALYSIS/alignment_filt/$STATION

	# Extract header
	samtools view -H $ANALYSIS/alignment/$STATION/$SAMPLE'.sort.bam' > $ANALYSIS/alignment_filt/$STATION/header.sam

	### Identidy filter

	# Create md bam
	samtools fillmd -b -e $ANALYSIS/alignment/$STATION/$SAMPLE'.sort.bam' $REF > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.md.bam'
	samtools view $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.md.bam' > $ANALYSIS/alignment_filt/$STATION/encodefile
	samtools view $ANALYSIS/alignment/$STATION/$SAMPLE'.sort.bam' > $ANALYSIS/alignment_filt/$STATION/completefile

	# 95% identity threshold
	$FILTER_ID $ANALYSIS/alignment_filt/$STATION/encodefile $ANALYSIS/alignment_filt/$STATION/completefile $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95'

	# Filtered bam
	cat $ANALYSIS/alignment_filt/$STATION/header.sam $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95' \
	  | samtools view -Sb - > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.bam'

	# Sort and index
	samtools sort $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.bam' > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.sort.bam'
	samtools index $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.sort.bam'

	### Low complexity filter

  # Create fasta file
	samtools bam2fq $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.sort.bam' | seqtk seq -A - > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.fasta'

	# Filter
	echo "PRINSEQ STEP "$STATION" SAMPLE "$SAMPLE
	prinseq-lite -fasta $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.fasta' -out_good $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc-good' -out_bad $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc-bad' -lc_method dust -lc_threshold 7 > $ANALYSIS/alignment_filt/'log_'$STATION'_'$SAMPLE'_low_complexity.txt' 2>&1

	samtools view $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_95.sort.bam' > $ANALYSIS/alignment_filt/$STATION/bam95

	grep ">" $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc-good.fasta' | sed 's/>/ @ /g' | sed 's/\// @ /g' | awk -F" @ " '{print $2}' | uniq > $ANALYSIS/alignment_filt/$STATION/good_read_lc
	grep ">" $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc-bad.fasta' | sed 's/>/ @ /g' | sed 's/\// @ /g' | awk -F" @ " '{print $2}' | uniq > $ANALYSIS/alignment_filt/$STATION/bad_read_lc

	NGOOD=$(wc -l $ANALYSIS/alignment_filt/$STATION/good_read_lc | cut -f1 -d' ')
	NBAD=$(wc -l $ANALYSIS/alignment_filt/$STATION/bad_read_lc | cut -f1 -d' ')


	if [ $NGOOD -le $NBAD ]
	then
		$FILTER_LCGOOD $ANALYSIS/alignment_filt/$STATION/bam95 $ANALYSIS/alignment_filt/$STATION/good_read_lc $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc95'
	else
		$FILTER_LCBAD $ANALYSIS/alignment_filt/$STATION/bam95 $ANALYSIS/alignment_filt/$STATION/bad_read_lc $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc95'
	fi
	echo "PRINSEQ FILTERING DONE "$STATION

  # New bam
	cat $ANALYSIS/alignment_filt/$STATION/header.sam $ANALYSIS/alignment_filt/$STATION/$SAMPLE'_lc95' \
	  | samtools view -Sb - > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filtered.bam'

	# Sort et Index
	samtools sort $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filtered.bam' > $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filter.sort.bam'
	samtools index $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filter.sort.bam'

	# Clean
	rm $ANALYSIS/alignment_filt/$STATION/encodefile $ANALYSIS/alignment_filt/$STATION/completefile $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.md.bam' $ANALYSIS/alignment/$STATION/$SAMPLE'.sort.bam' $ANALYSIS/alignment/$STATION/header.sam $ANALYSIS/alignment/$STATION/bad_read_lc $ANALYSIS/alignment/$STATION/bam95


	## Variant calling --------------
	mkdir $ANALYSIS/variants/$STATION

	# Coverages
	samtools mpileup -aa -B -A -f $ANALYSIS/variants/ref/$ASSEMBLY $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filter.sort.bam' > $ANALYSIS/variants/$STATION/$SAMPLE.mpileup

	tot=$(wc -l $ANALYSIS/variants/$STATION/$SAMPLE.mpileup | cut -f1 -d' ')
	hsum=$(cat $ANALYSIS/variants/$STATION/$SAMPLE.mpileup | awk '{if ($4 > 0) SUM += 1} END {print SUM}')
	vsum=$(cat $ANALYSIS/variants/$STATION/$SAMPLE.mpileup | awk '{SUM += $4} END {print SUM}')
	hcov=$(echo "scale=3; $hsum/$tot" | bc)
	vcov=$(echo "scale=3; $vsum/$tot" | bc)
	echo "$STATION $SAMPLE $tot $hsum $vsum $hcov $vcov" >> $ANALYSIS/coverage_$3

	# BCF file
	samtools mpileup --output-tags AD -D -uf $ANALYSIS/variants/ref/$ASSEMBLY $ANALYSIS/alignment_filt/$STATION/$SAMPLE'.filter.sort.bam' > $ANALYSIS/variants/$STATION/$SAMPLE'_rawcall.bcf'

	bcftools call --ploidy 1 -v -m $ANALYSIS/variants/$STATION/$SAMPLE'_rawcall.bcf' > $ANALYSIS/variants/$STATION/$SAMPLE'_calls.vcf'


	## Variant filtering--------------------------------------------------------------------------
	bcftools filter --exclude 'QUAL < 30' $ANALYSIS/variants/$STATION/$SAMPLE'_calls.vcf' | bcftools view -g ^miss > $ANALYSIS/variants/$STATION/$SAMPLE'_filtered_calls.vcf'
	bcftools view -v snps -O z $ANALYSIS/variants/$STATION/$SAMPLE'_filtered_calls.vcf' > $ANALYSIS/variants/$STATION/$SAMPLE'_snp.vcf.gz'
	bcftools index $ANALYSIS/variants/$STATION/$SAMPLE'_snp.vcf.gz'
done

# Fichier global
FILE=$(ls $ANALYSIS/variants/'TARA_'*/*'_snp.vcf.gz')
bcftools merge $FILE > $ANALYSIS/$3'_snp.vcf'
