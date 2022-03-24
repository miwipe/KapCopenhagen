# The Kap København (Cape Copenhagen) metagenomic analysis


This code are analyses that accompanies the Kjaer et al. 2022 article, and allows the reader to replicate the analysis in here. It consists of subdivisions of bash scripts and guides on how to:
1. Setting up proper environments using conda whenever possible. 
2. Download and prepare genomic databases to align DNA reads against. 
3. QC and mapping of all the reads against all databases.
4. Code for performing the taxonomic profiling, DNA damage and read length estimates, using metaDMG.
5. Rscript for parsing, filtering and plotting profiles.
6. Extraction of unique reads classified to focal taxa that was used for downstream analysis
7. PathPhynder placements of chloroplasts and mitochondrias
8. Plant phylogenetic placement and molecular dating
9. Mammalian phylogenetic placement and molecular dating
10. Performing the marine eukaryote SMAGs analysis


All code, database build, mapping, analysis and DNA damage estimates was performed on a Red Hat Enterprise Linux Server running 7.7 (Maipo), with 88 cores and 1Tb memory.


## Creating environment and installing dependencies using conda

### Making conda environment with provided environment.yml file
```
conda env create -f environment.yml
conda activate KapK
```
### Installing additional dependencies
```
mkdir programmes
cd programmes

wget http://kirill-kryukov.com/study/tools/fasta-splitter/files/fasta-splitter-0.2.6.zip 
unzip fasta-splitter-0.2.6.zip 

git clone https://github.com/keenerd/gz-sort; cd gz-sort; make; ./gz-sort -h

cd ../../
```

## Downloading and building the RefSeq databases


### Downloading fungi, viral and archaea reference genomes
```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*genomic.fna.gz
```
### Downloading other vertebrate references than mammalian

```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/*genomic.fna.gz
gzip -d *
cat *.fna > vert_other.fa
rm *.fna
/programmes/fasta-splitter.pl --n-parts 9 vert_other.fa
rm vert_other.fa
i=1
for file in vert_other.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in vert_other.?
do
bowtie2-build --threads 50 $file $file
done
```

### Downloading mammalian references 
``` 
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/*genomic.fna.gz
gzip -d *
cat *.fna > vert_mam.fa
rm *.fna
/programmes/fasta-splitter.pl --n-parts 18 vert_mam.fa
rm vert_mam.fa
i=1
for file in vert_mam.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in vert_mam.*
do
bowtie2-build --threads 50 $file $file
done
```

### Downloading invertebrate references 

```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/*genomic.fna.gz
gzip -d *
cat *.fna > invert.fa
rm *.fna
/programmes/fasta-splitter.pl --n-parts 3 invert.fa
rm invert.fa
i=1
for file in invert.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in invert.?
do
bowtie2-build --threads 50 $file $file
done
```

### Other available references, including mitochondria, plants, protozoans and plastids
```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz
gzip -d *
cat *.fna > others.fa
rm *.fna
/programmes/fasta-splitter.pl --n-parts 5 others.fa
rm others.fa
i=1
for file in others.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in others.?
do
bowtie2-build --threads 50 $file $file
done
```
### merging raw data per sample and trimming adaptors

Create a list of uniq sample names which needs to be merged cutting option -f might vary depending on your local file system. 
```
ll *fastq.gz | cut -f7 -d/ | uniq > merge.list

while read -r line
do
arr=($line)
lib=${arr[0]}
echo $lib
AdapterRemoval --file1 $lib*R1*.fastq.gz --file2 $lib*R2*.fastq.gz --mm 3 --minlength 30 --trimns --trimqualities --minquality 30 --basename $lib --collapse --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT 2> adap_rem.$lib.log.txt
done < merge.list

```

Next merge the collapsed and the collapsed truncated files per sample

```
while read -r line
do
arr=($line)
#if [ "${arr[1]}" = "$(basename $folder)" ]
lib=${arr[0]}
echo $lib
cat $lib.collapsed $lib.collapsed.truncated > $lib.col.fq &
done < merge.list
```

### QC filtering, mapping, merging and sorting alignments

All 

```
./holi.sh &> holi.log.txt
```



# Mammalian mitochondrial phylogenetic placement "Binia put you guide and/or code here, bash scripts in the scripts folder"


# Plant chloroplast phylogenetic placement "Bianca put you guide and/or code here, bash scripts in the scripts folder"


# Marine eukaryotic SMAGs analysis Antonios code




