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

### Making conda environment with the provided environment.yml file in the data folder
```
conda env create -f environment.yml
conda activate KapKBH
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
## Downloading and building the PhyloNorway arctic database
In your browser navigate to the PhyloNorway plant genome repo here https://dataverse.no/dataset.xhtml;jsessionid=dfe334bbadc5fb9c5eab4332d568?persistentId=doi:10.18710/3CVQAG&version=DRAFT
make sure you have enough data storage available as these file take up +200 GB storage, and will take up even more once index by bowtie2. 

Now in your terminal convert fasta file into indexed libraries, by:
``` 
for file in *fasta
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

# Molecular dating of the ancient Betula chloroplast from Kap København.

This workflow will describe how we performed the molecular dating analysis on the Betula chloroplast sequence.

### 1. Gather and align reference sequences.

Download reference sequences from NCBI. The samples used here were Betula chloroplast sequences:
```
KX703002.1 LC542973.1 MG386367.1 MG386368.1 MG386369.1 MG386370.1 MG386401.1 MG674393.1 MH205735.1 MK888853.1 MN830400.1 MT310900.1 MT872524.1 MT872525.1 MT872526.1 MT872527.1 MT872528.1 MT872529.1 MT872530.1 NC_033978.1 NC_037473.1 NC_039992.1 NC_039993.1 NC_039994.1 NC_039995.1 NC_039996.1 NC_047177.1
```
which includes a single Alnus sample, the closest outgroup with a tight chloroplast divergence date. This outgroup sample is necessary to calibrate the phylogeny.

We also used three extra sequences not on NCBI, from the PhyloNorway plant genome dataset (which can be found in the data folder). These were:
```
Betula_fruticosa-216991.TROM_V_201226.BXA-CDT.chloro.fasta
Betula_nana_nana-216990.TROM_V_96939.BMB-CC.chloro.fasta
Betula_nana_ssp_Tundrarum-216990.TROM_V_94022.BXA-CCC.chloro.fasta
```
for a total of 31 reference sequences.

Open these sequences in some annotation software (we used OGview, available online) to make sure they're all going in the same direction and start in the same place, otherwise they won't align well. Chloroplasts are circular, so it's not always the case that reference genomes will start at the same position and go in the same direction. We also made some pairwise dot plots (eg. using nucmer) to make sure they're all going in the same direction.

We wrote a homegrown program called "rotate" to make sure they all start at the same point and go the same direction, using an anchor string. It's available in the folder "code".

Align them with Mafft (or your program of choice) to get an MSA.
```
mafft --memsave betulas.fasta > betulas_mafft.afa
```
Open up this MSA in a viewer to check quality. We used NCBI's MSAviewer and Seqotron.

Call a consensus on the reference sequence MSA. Use Python BioAlign's consensus function. Python code file is get_consensus, available here.


### 2. Map ancient reads.
From the lca files in the metaDMG output (run in mode -simscoremin 0.95-100, and 0.90-0.95, we extracted readIDs classified to taxonomic level or lower:
```
for file in /data/lca/*lca.txt; do grep 'Betula' $file | awk -F":" '{print ">"$1":"$2":"$3":"$4":"$5":"$6":"$7 }' > $file.Betula.readID.txt & done
```
We then used the subseq option in seqtk to extract all read info from each individual fastq file (below $fq) that matched the readIDs extracted from the corresponding lca file:
```
seqtk subseq $file $file.Betula.readID.txt  > $file.Betula.readID.fq
```
We next mapped the fastqs with classified Betula reads using BWA with ancient DNA settings against a consensus Betula chloroplast, generated by aligning all chloroplast from the genus Betula and setting all variants as N's ## Bianca you might have the command for that?
```
for file in *Betula*fq
do
bwa aln -l 1024 -o 2 -n 0.001 -t 20 /willerslev/edna/KapKBH/cpDNA_Phylogenies/Consensus_dbs/betula_consensus.fa $file > $file.sai
bwa samse /willerslev/edna/KapKBH/cpDNA_Phylogenies/Consensus_dbs/betula_consensus.fa $file.sai $file > $file.sam
done
```
Remove unmapped reads, quality filter for read quality 25, sort:
```
for file in *; do name="${file%".fq.sam"}"; samtools view -F 0x04 -b $file > $name.aligned.bam; samtools view -q 25 -b $name.aligned.bam > $name.q25.bam; samtools sort $name.q25.bam > $name.q25.sorted.bam; done
```
End up with:
```
119_B3_116.Betula.q25.sorted.bam 
50_B3_127.Betula.q25.sorted.bam 
69_B2_100.Betula.q25.sorted.bam 
69_B2_103.Betula.q25.sorted.bam 
69_B2_97.Betula.q25.sorted.bam 
74_B1_83.Betula.q25.sorted.bam 
75_B1_83.Betula.q25.sorted.bam
```

In our case, it is appropriate to consider these as one sample. Merge these with samtools merge, and use samtools to sort, to get MergedBetula.q25.sorted.bam. Thus increasing breadth of coverage and depth of coverage.

Now we need to turn this into a fasta file, because that's what Beast takes as input. Use:
```
bcftools mpileup -f betula_consensus.fa MergedBetula.q25.sorted.bam > MergedBetula.q25.sorted.mpileup
bcftools call --ploidy 1 -m -P 0 -Ov MergedBetula.q25.sorted.mpileup > MergedBetula.q25.vcf
```
Save the mpileup so that we can pull out the per-site depths later. This will be an important filtering criterion to make sure we only use reliable sites.

The call options disable the default calling algorithm, which biases SNP calls towards a prior (in this case, the consensus sequence). Instead, we force a strict majority call and avoid this bias.
Notably, the default mpileup option uses a base alignment quality algorithm and then a base quality cutoff of 13. We experimented with turning this off, but found it to be beneficial, in that it greatly reduced false SNP calls near indel regions.

The following code makes the ancient vcf into a fasta file by replacing the sites in the consensus with the ancient genotypes, where available. It does so only for SNVs, because we don't believe indels in our ancient sequence to be sufficiently reliable to use in this kind of analysis.
```
bgzip -c MergedBetula.q25.vcf > MergedBetula.q25.vcf.gz
tabix -p vcf MergedBetula.q25.vcf.gz
bcftools view --types snps MergedBetula.q25.vcf.gz > MergedBetula.q25.snps.vcf
bgzip -c MergedBetula.q25.snps.vcf > MergedBetula.q25.snps.vcf.gz
tabix -p vcf MergedBetula.q25.snps.vcf.gz
bcftools consensus --missing N --sample MergedBetula.q25.sorted.bam -f betula_consensus.fa MergedBetula.q25.snps.vcf.gz > MergedBetula.snps.fa
```
Since the ancient sequence doesn't cover every region of the consensus sequence (in our case, some regions have no coverage at all), this fasta file will of course contain regions belonging to the consensus and not the ancient sequence. We will cut these regions out later.
```
cat betulas_mafft.afa MergedBetula.snps.fa > all_betulas.afa
```

### 3. Quality filtering, sanity checks, and making partitions for Beast.

Much of this was done manually.

Different types of sites evolve at different rates. Get the gff3 annotation file for the longest reference sequence, MG386368, from NCBI. We used the gff3 and fasta files for MG386368 to label indivudal sites as (1) coding regions, in which case we label as the position in the codon (1, 2 or 3, dependent on both the phase and strand noted in the gff3 file), (2) RNA, (3) neither coding nor RNA. We also pulled out the coding sequences into a separate file, reverse complemented and shifted where necessary according to the noted phase and strand, with each CDS region split by a "---". We want to double check the correctness of this annotation by making sure that these genes translate well, eg. with no premature stop codons, and that our labelling of the codon positions is correct. To do this, we opened up the coding sequence file in Seqotron (though any sequence viewer that performs translation can be used) and translated the sequence into proteins.

Since the ancient fasta file is on the same coordinates as the MSA, we could push the annotation through the multiple sequence alignment, then pull out the coding regions in the ancient sequence by using the annotation for MG386368. With no depth cutoff whatsoever, translating these regions leads to some premature stop codons and a poor quality protein alignment. However, when increasing our depth cutoff to 20 (by using the depths in the mpileup created above) and replacing sites in the ancient fasta which did not meet this filter by "N"s, we see a good protein alignment (except for the Ns, of course). This tells us that our ancient sample, at coding sites with at least 20 depth, are highly reliable. We experimented with different depth cutoffs here, and chose the lowest value that gave a reliable translation.

Another good thing to check: Where in the genome does the ancient fasta differ from the consensus, and from other sequences in the alignment? Look at both the protein sequence and the full sequence. Are there any regions where the ancient sequence seems to have many differences in a small region? If so, investigate these further. In our case, for example, we found 5 differences in a 100 base window. It turned out that the region was highly similar to a part of the Betula mitochondrion. This region was cut out with a depth cutoff of 20.

Next, push the annotation through the multiple sequence alignment. Check how many polymorphic sites, both in the reference and the ancient sequence, sit in each partition, and see if this makes sense to you. For example, you would expect a higher nucleotide substitution rate in the 3rd base of coding regions than in the first or second. This will inform how to make the partitions for Beast - in particular, we want to group similar regions with sufficient polymorphisms to be informative.

In our case, we used the following partitions. First, we only took sites which met a depth cutoff of at least 20 in the ancient sample, and which had less than three gaps "-"s in the multiple sequence alignment. Of these sites, we partitioned into (1) positions 1 or 2 in coding regions, (2) position 3 in coding regions, (3) RNA or (4) non-coding and non-RNA.

These regions contained 11668, 5828, 2690 and 29538 sites, respectively.


### 4. Run Beast.
Parameter choices here will also be highly dataset-dependent.

We used the four partitions described above, unlinked substitution models, a GTR+Gamma nucleotide substitution models with estimated base frequencies and 4 Gamma categories, and a strict clock. We assigned an age of 0 to the reference sequences, used a normal distribution prior with mean 61.1 million years and standard deviation 1.633 million years for the root height (see manuscript for reference; standard deviation was obtained by conservatively converting the 95% HPD to z-scores). We tried giving the ancient age several priors such as a Gamma (shape=1, scale=1.7) and a uniform, but found our results insensitive to the prior (see manuscript). Reported results are for a uniform prior with a minimum of 0 and a maximum of 10 million years. We also ran the coding regions alone, since they are able to translate correctly and therefore highly reliable sites, and found that they gave the same median and a much larger confidence interval, as expected when using fewer sites (see manuscript). We ran the MCMC for a total of 100 million iterations and used default options otherwise. We verified convergence in Tracer, and checked that the mutation rates of the different partitions were acting as expected. We also verified that the resulting MCC tree from TreeAnnotator had placed the ancient sequence phylogenetically in the same place as the pathPhynder placement.

We attach the multiple sequence alignment used and the main BEAST xml file.

This gave a 95% HPD of [-2.0172,-0.6786] and a median of 1.323 million years for the age of the ancient Betula chloroplast sequence.

# Marine eukaryotic SMAGs analysis Antonios code 
Link to repo? 




