### 1 Assembly MtDNA reads using CANU
```
canu \
   useGrid=false \
   -p Aamphitrite_mt \
   -d . \
   genomeSize=16k \
   -pacbio [filename] 
   # or -nanopore [filename]
```