### Run BUSCO with preffered lineages
```
TARGET=[genome_filename]
LABEL=[run_label]
CPU=[# cpus]

for db in [lineage_list]
do
  busco -i ${TARGET} -o ${LABEL}_${db} \
        -l ${db} -m genome -c $CPU \
        -f --offline --download_path [path_to_BUSCO_databases]
done
```