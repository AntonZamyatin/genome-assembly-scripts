### 1 Minimap overlap ONT reads over draft genome
```
CONTIGS=[draft_genome_filename]
READS=[ONT_reads_filename]
CPU=[# cpus]

minimap2 -x map-ont -t $CPU  $CONTIGS $READS > contigs_v_reads_ovl.paf
```
### 2 Polishing ONT assembly with ONT reads 1 round of Racon
```
racon -m 8 -x -6 -g -8 -w 500 -t $CPU $READS ./contigs_v_reads_ovl.paf $CONTIGS > contigs_racon.fasta
```