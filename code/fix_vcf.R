#!/usr/local/bin/Rscript


# May 20 2021, Bianca De Sanctis.
# This is a script to fix up the deletions/missing data represented by *s in the vcf files, eg. that you get from snp-sites output.
# Pathphynder can't cope with them. Also this will take out triallelic sites. 

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){print("Usage: Rscript fix_vcf.R input.vcf output.vcf. This code assumes your vcf has four lines in the header."); quit()}

vcfname = args[1]
outputname = args[2]

#############################
## should work for any vcf file ##

library(stringr)

# Get the header

header = strsplit(readLines(vcfname)[4],split="\t")[[1]]
original.vcf = read.table(vcfname,stringsAsFactors=FALSE,colClasses="character")
colnames(original.vcf) = header; numsites = nrow(original.vcf); ncols = ncol(original.vcf)
fixed.vcf = original.vcf 

lines.with.asterisks = grep("\\*",original.vcf$ALT)  # 358 lines have asterisks
print(paste0(length(lines.with.asterisks)," lines have asterisks, out of ",nrow(original.vcf)," total sites. Fixing..."))
lines.with.asterisks.first = grep("^\\*",original.vcf$ALT)
lines.with.asterisks.last = grep("\\*$",original.vcf$ALT)
# split into cases. you can look at all the possibilities with unique(original.vcf$ALT)
num.extra.alleles = str_count(original.vcf$ALT,",")

# case 1. biallelic (an alt and an asterisk)
biallelic.sites.with.asterisks = intersect(which(num.extra.alleles==1),lines.with.asterisks)
# case 1a. the asterisk is first.
biallelic.sites.with.asterisks.first = intersect(biallelic.sites.with.asterisks,lines.with.asterisks.first)
for(site in biallelic.sites.with.asterisks.first){
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==1)+9] = "."
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==2)+9] = "1"
  fixed.vcf[site,]$ALT = gsub("\\*,","",original.vcf[site,]$ALT)
}
# case 1b. the asterisk is second.
biallelic.sites.with.asterisks.second = setdiff(biallelic.sites.with.asterisks,lines.with.asterisks.first)
for(site in biallelic.sites.with.asterisks.second){
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==2)+9] = "."
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==1)+9] = "1" # no change, here for consistency
  fixed.vcf[site,]$ALT = gsub(",\\*","",original.vcf[site,]$ALT)
}
# case 2. triallelic (two alts and an asterisk)
triallelic.sites.with.asterisks = intersect(which(num.extra.alleles==2),lines.with.asterisks)
# case 2a. asterisk is first.
triallelic.sites.with.asterisks.first = intersect(triallelic.sites.with.asterisks,lines.with.asterisks.first)
for(site in triallelic.sites.with.asterisks.first){
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==1)+9] = "."
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==2)+9] = "1"
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==3)+9] = "2"
  fixed.vcf[site,]$ALT = gsub("\\*,","",original.vcf[site,]$ALT)
}
triallelic.sites.with.asterisks.second = setdiff(setdiff(triallelic.sites.with.asterisks,lines.with.asterisks.last),lines.with.asterisks.first)
for(site in triallelic.sites.with.asterisks.second){
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==1)+9] = "1" # consistency
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==2)+9] = "."
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==3)+9] = "2"
  fixed.vcf[site,]$ALT = gsub("\\*,","",original.vcf[site,]$ALT)
}
triallelic.sites.with.asterisks.third = intersect(triallelic.sites.with.asterisks,lines.with.asterisks.last)
for(site in triallelic.sites.with.asterisks.third){
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==1)+9] = "1" # consistency
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==2)+9] = "2" # consistency
  fixed.vcf[site,which(fixed.vcf[site,10:ncols]==3)+9] = "."
  fixed.vcf[site,]$ALT = gsub(",\\*$","",original.vcf[site,]$ALT)
}

# ok good, we've gotten rid of all the *s. now let's look at remaining triallelic sites. 
# If one of them is a G->A or C->T this is probable deamination and we want to mark it missing, especially in cases where there are not very many of them with the deaminated allele.
# i haven't automated this part though. as of now i'm just deleting them. 
# if you want to cope with this in some way that isn't just deleting the sites,
# it will need different code.
not_biallelic = which(str_count(fixed.vcf$ALT,",")>=1)
print(paste0(length(not_biallelic)," sites are not biallelic, even after fixing asterisks. These will be deleted."))
# too bad... i guess there's a lot of variation... there are 433 triallelic sites.
# but pathphynder just can't cope with those. 
# so i guess we will take them out.
fixed.vcf = fixed.vcf[-not_biallelic,]

print(paste0("The fixed vcf has ",nrow(fixed.vcf)," biallelic sites.")) # this is how many snps you're left with (biallelic ones)

# write the fixed vcf to file
conn = file(outputname,'w')
writeLines(readLines(vcfname)[1:4],con=conn)
write.table(fixed.vcf,file=conn,append=TRUE,quote=FALSE,
            row.names=FALSE,col.names=FALSE,sep="\t")








