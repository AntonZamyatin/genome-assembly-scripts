### 1 round Racon + Medaka polishing ONT assembly with ONT reads
```
CONTIGS=[draft_genome_filename]
READS=[ONT_reads_filename]
CPU=[# cpus]

minimap2 -x map-ont -t $CPU  $CONTIGS $READS > contigs_v_reads_ovl.paf

racon -m 8 -x -6 -g -8 -w 500 -t $CPU $READS ./contigs_v_reads_ovl.paf $CONTIGS > contigs_racon.fasta

CONTIGS=contigs_racon.fasta

minimap2 -ax map-ont -t $CPU $CONTIGS $READS > tmp.sam

samtools view -@ $CPU -bhS tmp.sam > aln.bam
samtools sort -@ $CPU aln.bam -o aln_sorted.bam
samtools index aln_sorted.bam
rm rmp.sam
rm aln.bam

PARAMS=r941_prom_high_g4011 # pore and basecaller params

mkdir medaka_consensus
medaka_consensus -m $PARAMS -i $READS -d $CONTIGS -o medaka_consensus -t $CPU
```