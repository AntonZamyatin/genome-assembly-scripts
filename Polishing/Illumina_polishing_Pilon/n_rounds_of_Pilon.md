### Full cycle of n rounds of Pilon polishing
```
READS1=[fwd_reads_filename]
READS2=[rvs_reads_filename]
CURDRAFT=[draft_genome_filename]
CPU=[# cpus]

for n in 1 2 3
do
    WORKDIR=pilon${n}
    mkdir $WORKDIR
    cp $CURDRAFT ${WORKDIR}/draft${n}.fasta
    TARGET=${WORKDIR}/draft${n}.fasta

    bwa index ${TARGET}

    bwa mem ${TARGET} ${READS1} ${READS2} -t $CPU | \
        samtools view -bS > ${WORKDIR}/aln.bam

    samtools sort ${WORKDIR}/aln.bam -o ${WORKDIR}/aln_sorted.bam
    samtools index ${WORKDIR}/aln_sorted.bam
    rm ${WORKDIR}/aln.bam

    pilon -Xmx160048m --genome ${TARGET} --bam ${WORKDIR}/aln_sorted.bam --outdir ${WORKDIR} \
          --diploid --threads $CPU \
          --changes --fix bases > ${WORKDIR}/pilon${n}.out

    python pilonOutPatser -i ${WORKDIR}/pilon${n}.out -o ${WORKDIR}/pilon${n}.stats

    CURDRAFT=${WORKDIR}/pilon.fasta
done
```
