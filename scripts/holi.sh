
for infile in $(pwd)/*.fq
do
bname=$(basename $infile)
echo $bname
bname2=$(echo $bname | sed 's/.fq*/_holi/')
basepath=$(pwd)/
basefolder=$basepath
echo $basepath
echo $bname2
mkdir $basepath$bname2
cd $basepath$bname2
pwd

echo Step 1. Removing poly A tails
fastq-grep -v "AAAAA$" ../$bname > kmer_$bname
echo Step 2. Removing reverse complemented A tails
fastq-grep -v "^TTTTT" kmer_$bname > kmer2_$bname
echo Step 3. Removing rememnants adapter sequence 1 = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" kmer2_$bname > adap1_kmer2_$bname
echo Step 4. Removing remnants adapter sequence 2 = ATCTCGTATGCCGTCTTCTGCTTG
fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" adap1_kmer2_$bname > adap2_kmer2_$bname


echo Step 6. sga preprocessing
sga preprocess --dust-threshold=1 -m 30 adap2_*.fq -o adap2_kmer2_$bname.pp.fq
echo Step 7. sga index
sga index --algorithm=ropebwt --threads=30 adap2_kmer2_$bname.pp.fq
echo Step 8. sga filter
sga filter --threads=30  --no-kmer-check adap2_kmer2_$bname.pp.fq -o adap2_kmer2_$bname.pp.rmdup.fq
echo Step 9. Calculating read length distribution and outputting file
cat adap2_kmer2_$bname.pp.rmdup.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > adap2_kmer2_$bname.pp.rmdup.fq.read_length.txt

ls -lh
rm kmer*
rm adap1*
rm adap2_kmer2_$bname
rm adap2_kmer2_$bname.pp.fq
rm *wt
rm *sai
rm *rmdup.discard.fa

for DB in /database/norPlantCom.?
do
echo Mapping adap2_kmer2_$bname.pp.rmdup.fq against $DB
nice -20 bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/nt.?
do
echo Mapping adap2_kmer2_$bname.pp.rmdup.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/vert_other.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/vert_mam.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/vert_mam.??
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/invert.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/viral_fungi_archaea.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/gtdb_r89
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/arctic_animals.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/arctic_animals_other.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /database/arctic_animal_b3.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done


samtools merge --verbosity 5  $bname.merged.sam.gz *.bam -@ 60

ls -lh *bam

rm *nt.?*
rm *vert_other.?*
rm *vert_mam.?*
rm *vert_mam.??*
rm *invert.?*
rm *norPlantCom*
rm *viral_fungi_archaea*
rm *arctic_animal*

cd $basepath
done
