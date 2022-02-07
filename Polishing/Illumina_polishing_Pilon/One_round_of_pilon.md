### 1 BWA alignment Illumina reads over draft genome
```
READS1=[fwd_reads_filename]
READS2=[rvs_reads_filename]
# ILL_FOLDER=/lustre/groups/cbi/Projects/barnacles/Aamphitrite/data/genomic/ILLUMINA_DATA/
WORKDIR=pilon_out
TARGET=[draft_genome_filename]
CPU=[# cpus]

bwa index ${TARGET}

bwa mem ${TARGET} ${READS1} ${READS2} -t $CPU | \
        samtools view -bS > ${WORKDIR}/aln.bam
```
### 2 Sorting and indexing alignment
```
samtools sort ${WORKDIR}/aln.bam -o ${WORKDIR}/aln_sorted.bam
samtools index ${WORKDIR}/aln_sorted.bam
rm ${WORKDIR}/aln.bam
```
### 3 Pilon polishing
```
pilon -Xmx160048m --genome ${TARGET} --bam ${WORKDIR}/aln_sorted.bam --outdir ${WORKDIR} \
    --diploid --threads 40 \
    --changes --fix bases > pilon.out
```
### 4 Parse pilon stats from outfile (parser script included in this folder)
```
python pilonOutPatser -i $WORKDIR/pilon.out -o $WORKDIR/pilon.stats
```
