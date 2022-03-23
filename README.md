# KapCopenhagen

This code are analyses that accompanies the Kjaer et al. 2022 article, and allows the reader to replicate the analysis in here. It consists of subdivisions of bash scripts and guides on how to:
1. Setting up proper environments using conda whenever possible. 
2. Download and prepare genomic databases to align DNA reads against. 
3. QC and mapping of all the reads against all databases.
4. Code for performing the taxonomic profiling, DNA damage and read length estimates, using metaDMG.
5. Rscript for parsing, filtering and plotting profiles.
6. Extraction of unique reads classified to focal taxa that was used for downstream analysis
7. Plant phylogenetic placement and molecular dating
8. Mammalian phylogenetic placement and molecular dating
9. Marine eukaryote SMAGs analysis


## Creating environment and installing dependencies using conda

```
conda env create -f environment.yml
conda activate KapK
```
```
mkdir programmes
cd programmes

wget http://kirill-kryukov.com/study/tools/fasta-splitter/files/fasta-splitter-0.2.6.zip 
unzip fasta-splitter-0.2.6.zip 

git clone https://github.com/keenerd/gz-sort; cd gz-sort; make; ./gz-sort -h

cd ../../
```

## Downloading and building databases

```

```


Database build, mapping and DNA damage estimates was performed on a Red Hat Enterprise Linux Server running 7.7 (Maipo). 
