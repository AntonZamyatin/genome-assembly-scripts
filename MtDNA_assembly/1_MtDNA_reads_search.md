### 1 Reads alignment to reference MtDNA sequence
```
REFERENCE=[filename]
READS=[filename]

# if pacbio reads:
ALN_PARAM=map-pb
# if ONT reads:
ALN_PARAM=map-ont
minimap2 -ax $ALN_PARAM $REFERENCE $READS > aln.sam
```

### 2 Filtration of mapped reads:
```
samtools view -F 4 aln.sam -@ 40 > mapped.sam
```

### 3 Get mapped reads from alignment .sam:
```
samtools fastq  mapped.sam > mapped_reads.fastq
```
