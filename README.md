# *Pectobacterium* and microbiome analysis of potato stem samples
Please download or obtain access to the following before doing the analysis

##  Programs
Clean metagenomic raw reads
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)

Metagenomic analysis 
- [Kraken2](http://ccb.jhu.edu/software/kraken/)
- [Sourmash](https://sourmash.readthedocs.io/en/latest/)

Map host genome and remove mapped reads
- [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
- [Samtools](https://www.htslib.org/doc/samtools.html)

## Illumina QC
### Trimmomatic
Trimmomatic cleans the metagenomic Illumina reads based on quality control. Trimmomatic uses Java, therefore for our script, we use the module system to use it. 
```bash
module load jdk/1.8.0
for i in $(ls *.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); java -jar trimmomatic-0.39.jar PE -phred33 \
#Seperate files into R1 and R2
${i}_R1.fastq ${i}_R2.fastq \
/cleanreads/$N.R1.paired.fastq
/cleanreads/$N.R1.unpaired.fastq \
/cleanreads/$N.R2.paired.fastq
/cleanreads/$N.R2.unpaired.fastq \
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100;done

#We use /cleanreads/ as the outputfolder
#make sure your PATH for your raw data does not contain underscores as that will be cut and the script will not work.
```

## Metagenome Analysis
### Kraken - Create a database using NCBI genomes
1. Download Taxonomy to create a database called DB
```bash
cd /kraken2/
./kraken2-build --download-taxonomy --db /folder/DB
```
2. Add genomes to the database library
```bash
#keep all NCBI genomes in the same folder (genomes/) to create a library with a loop.

for files in /genomes/*.fna; do \
./kraken2-build --add-to-library $file --db /folder/DB; done
```
3. Build the database
```bash
./kraken2-build --build --db /folder/DB
```

### Kraken - run the analysis
```bash
#input-location/ is a folder that contains all the clean metagenomic raw reads.

cd /kraken2
for i in $(ls /input-location/*.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); ./kraken2  \
--db /folder/DB \
--paired --report $N.report \
--output $N.kraken \
${i}_R1.fastq ${i}_R2.fastq; done
```

### Sourmash
Create a collection of signatures (a database)
1. Sourmash sketch to create the signature
```bash
#genomes is a folder that 
sourmash sketch dna -p k=31 genomes/*.fna -outdir genomesDB
```
2. Index the DB
```bash
sourmash index -k 31 genomesDB.sbt.zip genomesDB/*.sig
```
Create signatures for the metagenomic reads
```bash
for i in $(ls metagenomic_reads/*.fq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fq);sourmash sketch dna -p k=31 ${i}_R1.fq ${i}_R2.fq --name $N \
-o $N.zip;done
```
Compare DB signatures with the metagenomics reads
```bash
for i in $ $N.zip; do N=$(basename $i .zip); sourmash gather -p k=31 $i genomesDB.sbt.zip
```
## Remove host reads before metagenomic analysis
Another option is to remove the host (potato) reads before the analysis

### BWA-mem
```bash
module load bwa/0.7.17
bwa index hostgenome.fa \

for i in $(ls metagenomicreads/*.fastq | grep "_R1" | cut -f 1 -d "_"); do \
N=$(basename $i .fastq); bwa mem hostgenome.fa ${i}_R1.fastq ${i}_R2.fastq \
> $N.sam; done
```
### Samtools
The unmapped reads to the potato genome can be extracted from the SAM file. The unmapped reads correspond to the microbial reads.
1. convert to bam files and sort
```bash
cd /output-location
for i in /bwamem-output-location/*.sam; do \
N=$(basename $i .sam); samtools sort $i -O bam -o $N.bam; done
 ```   
2. Keeping unmapped reads
```bash
for i in *.bam; do \
N=$(basename $i .bam); samtools view -b -f 4 $i -o $N.unmapped.bam; done
 ```
3. sort unmapped
```bash
for i in *.unmapped.bam; do \
N=$(basename $i .unmapped.bam); samtools sort $i -n -o $N.sorted.bam; done
 ```
4. convert bam to fastq files
```bash
for i in *.sorted.bam; do \
N=$(basename $i .sorted.bam); samtools fastq $i > $N.nonpotato.fastq; done
```
