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
module load jdk/1.8.0
for i in $(ls /input-location/*.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); java -jar /trimmomatic-location/trimmomatic-0.39.jar PE -phred33 \
#Seperate files into R1 and R2
${i}_R1.fastq ${i}_R2.fastq \
/output-location/$N.R1.paired.fastq
/output-location/$N.R1.unpaired.fastq \
/output-location/$N.R2.paired.fastq
/output-location/$N.R2.unpaired.fastq \
ILLUMINACLIP:/trimmomatic-location/adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100;done

#change 'output/input/etc-location' to where you want to store the database on your device

```

## Metagenome Analysis
### Kraken - Create a database using NCBI genomes
```bash
#change 'output/input/etc-location' to where you want to store the database on your device

cd /kraken2
./kraken2-build --download-taxonomy --db /output-location

#add genomes to database library
for files in /genome-location/*.fna; do \
./kraken2-build --add-to-library $file --db /output-location; done

#build the database
./kraken2-build --build --db /output-location
```

### Kraken - run the analysis
```bash
#change 'output/input/etc-location' to where you want to store the database on your device

cd /kraken
for i in $(ls /input-location/*.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); ./kraken2  \
--db /kraken2-location/k2_standard_20230605 \
--paired --report /output-location/$N.report \
--output /output-location/$N.kraken \
${i}_R1.fastq ${i}_R2.fastq; done
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
module load bwa/0.7.17
bwa index /host-reference-genome-location/hostgenome.fa \
cd /host-reference-genome-location
for i in $(ls /input-location/*.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); bwa mem hostgenome.fa ${i}_R1.fastq ${i}_R2.fastq \
> /output-location/$N.sam; done
```
### Samtools
```bash
code
```
