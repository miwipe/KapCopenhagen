# The Kap København (Cape Copenhagen)


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



# Mammalian mitochondrial phylogenetic placement
## Concatenating and aligning mitogenome reference sequences downloaded from NCBI (for PathPhynder)
```
cat *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references.fa

mafft --thread n cat_NCBI_mitogenome_references.fa > Aln_NCBI_mitogenome_references.fa
```

## Create consensus sequence for sample (for BEAST)
```
angsd -dofasta 2 -docounts 1 -minmapq 25 -minq 25 -uniqueonly 1 -mininddepth 5 -i Sample.taxa.sort.bam -out Cons_Sample.taxa.depth5

gunzip Cons_Sample.taxa.depth5.fa.gz
bash rename_fasta_header.sh
```
## Concatenating and aligning all mitogenome reference sequences downloaded from NCBI and sample consensus (for BEAST)
```
cat Cons_Sample.taxa.depth5.fa *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references_query.fa
mafft --thread n cat_NCBI_mitogenome_references_query.fa > Aln_NCBI_mitogenome_references_query.fa
```
## Build consensus sequence for mitogenome reference sequences downloaded from NCBI
1.) Aln-file opened in Geneious, consensus sequences created with 75% Majority rule for family level/each clade

2.) Alignment created from all clade-consensus mitogenome references in Geneious

## Running BEAST (Phylogenetic placement mtDNA)
1.) Aln_NCBI_mitogenome_references_query.fa opened in BEAUti (v1.10.4)

2.) Run with default parameters, chain 2mio
```
beast2 -threads n Aln_NCBI_mitogenome_references_query.xml
```
3.) Tracer (v1.7.2) checked for convergence

4.) TreeAnnotator (v1.10.4), 10% Burnin removed, Maximum Clade Credibility Tree created

5.) FigTree (v1.4.4), Tree visualized with posterior probabilities

## Running BEAST (Reference tree for PathPhynder)
1.) Aln_NCBI_mitogenome_references.fa opened in BEAUti (v1.10.4)

2.) Run with default parameters, chain 2mio
```
beast2 -threads n Aln_NCBI_mitogenome_references.xml
```
3.) Tracer (v1.7.2) checked for convergence

4.) TreeAnnotator (v1.10.4), 10% Burnin removed, Maximum Clade Credibility Tree created

5.) FigTree (v1.4.4), Tree converted into Newick format for PathPhynder: ref_tree.nwk

## For PathPhynder (Phylogenetic placement mtDNA)
Input files:
- ref_tree.nwk
- Aln_NCBI_mitogenome_references.fa
- Query: sequencing_reads.fastq, e.g. Sample.taxa.fq

### For PathPhynder: Create VCF with snp-sites from multiple sequence alignment (MSA) file
```
snp-sites -v -c -o NCBI_mitogenome_references.vcf NCBI_mitogenome_references.fa
```
### For PathPhynder: Change missing data coding in vcf with Rscript
```
Rscript fix_vcf.R NCBI_mitogenome_references.vcf NCBI_mitogenome_references_fixed.vcf
```
### For PathPhynder: Fix consensus naming problem in vcf file (replace 1’s with “consensus” in first column of vcf) ###
```
awk '{ if ($1 == "1") $1="consensus";}1' NCBI_mitogenome_references_fixed.vcf | sed 's/ /\t/g' > NCBI_mitogenome_references_fixed_cons.vcf
```
### For PathPhynder: Create consensus for MSA.fa and index it ###
```
python get_consensus NCBI_mitogenome_references.fa Cons_NCBI_mitogenome_references.fa
bwa index Cons_NCBI_mitogenome_references.fa
```
### For PathPhynder: Mapping sample reads to reference consensus sequence ###
Create directory outputfolderpath/pathPhynder_analysis/map_to_cons/
```
bwa aln -l 1024 -n 0.001 -t 10 /path/to/Cons_NCBI_mitogenome_references.fa /path/to/samples/from/ngsLCA/Sample.taxa.fq | bwa samse /path/to/Cons_NCBI_mitogenome_references.fa  - /path/to/samples/from/ngsLCA/Sample.taxa.fq | samtools view -F 4 -q 25 -@ 10 -uS - | samtools sort -@ 10 -o outputfolderpath/path/Phynder_analysis/map_to_cons/Sample.taxa.sort.bam
```
Count number of reads mapped
```
samtools view -c file.sort.bam
```
Coverage estimates
```
bamcov -H file.sort.bam
```
Read depth across mitogenome
```
samtools depth file.sort.bam > file.sort.coverage
```
Read length distribution for mapped reads
```
samtools view file.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > file.readlength_dist.txt
```
### Running PathPhynder
Create directory pathphynder_results in e.g. outputfolderpath/pathPhynder_analysis/

### Running PathPhynder: 1.) Assigning informative SNPs to tree branches
```
/home/kbt252/Software/phynder/phynder -B -o /path/to/pathphynder_results/branches.snp /path/to/tree/ref_tree.nwk /path/to/fixed_cons_vcf/NCBI_mitogenome_references_fixed_cons.vcf
```
### Running PathPhynder: 2.) Call SNPs in a given dataset of ancient samples and find the best path and branch where these can be mapped in the tree
Prepare data: this will output a bed file for calling variants and tables for pylogenetic placement. Start command in results folder (e.g. pathphynder_results), as this command creates a new directory called “tree_data” in the folder you run the command in
```
pathPhynder -s prepare -i /path/to/tree/ref_tree.nwk -p taxa_pathphynder_tree -f /path/to/pathphynder_results/branches.snp  -r /path/to/Cons_NCBI_mitogenome_references.fa
```
### Running PathPhynder: 3.) Find best branch path and create tree file

FOR SINGLE BAM FILE TRANSITIONS AND TRANSVERSIONS
```
pathPhynder -s all -t 100 -i /path/to/tree/ref_tree.nwk -p /path/to/pathphynder_results/tree_data/taxa_pathphynder_tree -b outputfolderpath/pathPhynder_analysis/map_to_cons/Sample.taxa.sort.bam -r /path/to/Cons_NCBI_mitogenome_references.fa
```

FOR BAMLIST
```
pathPhynder -s all -t 100 -i /path/to/tree/ref_tree.nwk -p /path/to/pathphynder_results/tree_data/taxa_pathphynder_tree -l outputfolderpath/pathPhynder_analysis/map_to_cons/bamlist.txt -r /path/to/Cons_NCBI_mitogenome_references.fa
```

FOR TRANSVERSIONS ONLY
```
pathPhynder -s all -t 100 -m transversions -i /path/to/tree/ref_tree.nwk -p /path/to/pathphynder_results/tree_data/taxa_pathphynder_tree -b outputfolderpath/pathPhynder_analysis/map_to_cons/Sample.taxa.sort.bam -r /path/to/Cons_NCBI_mitogenome_references.fa
```

# Plant chloroplast phylogenetic placement "Bianca put you guide and/or code here, bash scripts in the scripts folder"


# Marine eukaryotic SMAGs analysis Antonios code




