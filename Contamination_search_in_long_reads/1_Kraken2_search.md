### Search contamination in reads over Kraken2 database
```
kraken2 --threads [num_cpu] -db [DB_path] [reads] \
    --classified-out calssified.out \
    --unclassified-out unclassified.out
```