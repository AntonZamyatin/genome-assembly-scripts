### Quast-lg run
```
TARGET=[genome_filename]
LABEL=[run_label]
OUTPUT="./${LABEL}"
REFERENCE=[reference_genome_filename]
CPU=40

mkdir $OUTPUT

quast.py -o ${OUTPUT} --threads ${CPU}\
         --labels ${LABEL} \
         --large -k \
         -r {$REFERENCE} \
         --eukaryote \
         ${TARGET}
```