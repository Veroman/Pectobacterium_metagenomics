# Pectobacterium and microbiome analysis of potato stem samples
Please download or obtain access to following before doing the analysis

##  Programs
Clean raw reads
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)

Metagenomic analysis 
- [Kraken](http://ccb.jhu.edu/software/kraken/)
- [Sourmash](https://sourmash.readthedocs.io/en/latest/)

Map host genome and remove mapped reads
- [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
- [Samtools](https://www.htslib.org/doc/samtools.html)

## Illumina QC
### Trimmomatic
```bash
#Trimmomatic is for quality control of Illumina reads
code
```

## Metagenome Analysis
### Kraken - Create a database using NCBI genomes
```bash
code
```

### Kraken - run the analysis
```bash
code
```
### Sourmash
Sourmask sketch to DB
```bash
sourmash sketch dna -p k=31 genomes/*.fna
```
Index to DB
```bash
sourmash index -k 31 CBDB signatures/*.sig
```
Sourmash sketch to metagenome
```bash
sourmash sketch dna -p k=31 contigs.fasta
```
Compare
```bash
sourmash gather contigs.fasta.sig 
```
## Remove host reads before metagenomic analysis
### BWA-mem
```bash
code
```
### Samtools
```bash
code
```
