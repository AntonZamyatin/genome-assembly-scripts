### 1 Minimap overlap ONT reads over draft genome
```
CONTIGS=contigs_racon.fasta
READS=haplotype-ZANU.fasta
CPU=[# cpus]

minimap2 -ax map-ont -t $CPU $CONTIGS $READS > tmp.sam
```

### 2 Alignment sorting and indexing
```
samtools view -@ $CPU -bhS tmp.sam > aln.bam
samtools sort -@ $CPU aln.bam -o aln_sorted.bam
samtools index aln_sorted.bam
rm rmp.sam
rm aln.bam
```

### 3 Medaka polishing
```
PARAMS=r941_prom_high_g4011 # pore and basecaller params

mkdir medaka_consensus
medaka_consensus -m $PARAMS -i $READS -d $CONTIGS -o medaka_consensus -t $CPU
```
