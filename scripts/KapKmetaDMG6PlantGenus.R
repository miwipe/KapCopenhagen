library(tidyverse) 
library(gplots)
library(scales)
library(reshape2)
library(vegan)
library(rioja)
library(readxl)
library(ggridges)
library(viridis)
library(hrbrthemes)

# import data
df <- read_csv("KapKmetaDMG6out.csv")

## changing all sample names
df$sample <- str_replace(df$sample, "KapK-12-1-24-Ext-1-Lib-1-Index2", "119_B3_116_L0_KapK_12_1_24")
df$sample <- str_replace(df$sample, "KapK-12-1-25-Ext-2-Lib-2-Index1", "119_B3_116_L0_KapK_12_1_25")
df$sample <- str_replace(df$sample, "KapK-12-1-27-Ext-4-Lib-4-Index2", "119_B3_116_L0_KapK_12_1_27")
df$sample <- str_replace(df$sample, "KapK-12-1-29-Ext-6-Lib-6-Index2", "119_B3_116_L0_KapK_12_1_29")
df$sample <- str_replace(df$sample, "KapK-12-1-31-Ext-8-Lib-8-Index2", "69_B2_97_L0_KapK_12_1_31")
df$sample <- str_replace(df$sample, "KapK-12-1-33-Ext-10-Lib-10-Index2", "69_B2_97_L0_KapK_12_1_33")
df$sample <- str_replace(df$sample, "KapK-12-1-34-Ext-11-Lib-11-Index2", "69_B2_100_L0_KapK_12_1_34")
df$sample <- str_replace(df$sample, "KapK-12-1-35-Ext-12-Lib-12-Index2", "69_B2_100_L0_KapK_12_1_35")
df$sample <- str_replace(df$sample, "KapK-12-1-36-Ext-13-Lib-13-Index2", "69_B2_100_L0_KapK_12_1_36")
df$sample <- str_replace(df$sample, "KapK-12-1-37-Ext-17-Lib-17-Index1", "69_B2_103_L0_KapK_12_1_37")
df$sample <- str_replace(df$sample, "KapK-12-1-38-Ext-18-Lib-18-Index1", "69_B2_103_L0_KapK_12_1_38")
df$sample <- str_replace(df$sample, "KapK-12-1-39-Ext-19-Lib-19-Index1", "69_B2_103_L0_KapK_12_1_39")
df$sample <- str_replace(df$sample, "KapK-12-1-41-Ext-20-Lib-20-Index1", "119_B3_116_L1_KapK_12_1_41")
df$sample <- str_replace(df$sample, "KapK-12-1-42-Ext-21-Lib-21-Index1", "119_B3_116_L1_KapK_12_1_42")
df$sample <- str_replace(df$sample, "KapK-12-1-43-Ext-22-Lib-22-Index1", "119_B3_116_L1_KapK_12_1_43")
df$sample <- str_replace(df$sample, "KapK-12-1-44-Ext-23-Lib-23-Index1", "119_B3_116_L2_KapK_12_1_44")
df$sample <- str_replace(df$sample, "KapK-12-1-45-Ext-24-Lib-24-Index1", "119_B3_116_L2_KapK_12_1_45")
df$sample <- str_replace(df$sample, "KapK-12-1-46-Ext-25-Lib-25-Index1", "119_B3_116_L2_KapK_12_1_46")
df$sample <- str_replace(df$sample, "KapK-12-1-47-Ext-26-Lib-26-Index1", "119_B3_116_L3_KapK_12_1_47")
df$sample <- str_replace(df$sample, "KapK-12-1-48-Ext-27-Lib-27-Index1", "119_B3_116_L3_KapK_12_1_48")
df$sample <- str_replace(df$sample, "KapK-12-1-49-Ext-28-Lib-28-Index1", "119_B3_116_L3_KapK_12_1_49")
df$sample <- str_replace(df$sample, "KapK-12-1-50-Ext-29-Lib-29-Index1", "119_B3_116_L4_KapK_12_1_50")
df$sample <- str_replace(df$sample, "KapK-12-1-51-Ext-33-Lib-33-Index2", "119_B3_116_L4_KapK_12_1_51")
df$sample <- str_replace(df$sample, "KapK-12-1-52-Ext-34-Lib-34-Index1", "119_B3_116_L4_KapK_12_1_52")
df$sample <- str_replace(df$sample, "KapK-198A-Ext-55-Lib-55-Index1", "50_B3_127_L0_KapK_198A")
df$sample <- str_replace(df$sample, "KapK-199A-Ext-56-Lib-56-Index1", "50_B3_127_L0_KapK_199A")
df$sample <- str_replace(df$sample, "KapK-199B-Ext-48-Lib-48-Index1", "50_B3_127_L0_KapK_199B")
df$sample <- str_replace(df$sample, "KapK-200A-Ext-49-Lib-49-Index1", "50_B3_127_L0_KapK_200A")
df$sample <- str_replace(df$sample, "KapK-202B-Ext-38-Lib-38-Index2", "74b_B1_83_L1_KapK_202B")
df$sample <- str_replace(df$sample, "KapK-203A-Ext-42-Lib-42-Index1", "74b_B1_83_L1_KapK_203A")
df$sample <- str_replace(df$sample, "KapK-203B-Ext-43-Lib-43-Index1", "74b_B1_83_L1_KapK_203B")
df$sample <- str_replace(df$sample, "KapK-205A-Ext-51-Lib-51-Index1", "74b_B1_83_L2_KapK_205A")
df$sample <- str_replace(df$sample, "KapK-205B-Ext-52-Lib-52-Index2", "74b_B1_83_L2_KapK_205B")
df$sample <- str_replace(df$sample, "KapK-205C-Ext-53-Lib-53-Index1", "74b_B1_83_L2_KapK_205C")
df$sample <- str_replace(df$sample, "KapK-205D-Ext-54-Lib-54-Index2", "74b_B1_83_L2_KapK_205D")
df$sample <- str_replace(df$sample, "Lok-75-Sample-1-Ext-58-Lib-58-Index1", "74a_B1_83_L1_Lok_75_Sample_1_Ext_58")
df$sample <- str_replace(df$sample, "Lok-75-Sample-1-Ext-A01-Lib01A-Index1", "74a_B1_83_L1_Lok_75_Sample_1_Ext_A01")
df$sample <- str_replace(df$sample, "Lok-75-Sample-1-Ext-A01-Lib01A-Index2", "74a_B1_83_L1_Lok_75_Sample_1_Ext_A01_2")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2a-Ext-A07-Lib11A-Index1", "74a_B1_83_L2_Lok_75_Sample_2a_Ext_A07")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2a-Ext-A08-Lib08A-Index2", "74a_B1_83_L2_Lok_75_Sample_2a_Ext_A08")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2a-Ext-A17-Lib17A-Index1", "74a_B1_83_L2_Lok_75_Sample_2a_Ext_A17")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2b-Ext-A09-Lib09A-Index1", "74a_B1_83_L2_Lok_75_Sample_2b_Ext_A09")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2b-Ext-A20-Lib20A-Index1", "74a_B1_83_L2_Lok_75_Sample_2b_Ext_A20")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2b-Ext-A21-Lib21A-Index1", "74a_B1_83_L2_Lok_75_Sample_2b_Ext_A21")
df$sample <- str_replace(df$sample, "Lok-75-Sample-2b-Ext-A22-Lib22A-Index2", "74a_B1_83_L2_Lok_75_Sample_2b_Ext_A22")
df$sample <- str_replace(df$sample, "Lok-75-Sample-3-Ext-59-Lib-59-Index2", "74a_B1_83_L3_Lok_75_Sample_3_Ext_59")
df$sample <- str_replace(df$sample, "Lok-75-Sample-3-Ext-A04-Lib04A-Index1", "74a_B1_83_L3_Lok_75_Sample_3_Ext_A04")
df$sample <- str_replace(df$sample, "Lok-75-Sample-3-Ext-A05-Lib05A-Index2", "74a_B1_83_L3_Lok_75_Sample_3_Ext_A05")
df$sample <- str_replace(df$sample, "Lok-75-Sample-3-Ext-A06-Lib06A-Index1", "74a_B1_83_L3_Lok_75_Sample_3_Ext_A06")
df$sample <- str_replace(df$sample, "Lok-75-Sample-4a-Ext-A11-Lib11A-Index1", "74a_B1_83_L4_Lok_75_Sample_4a_Ext_A11")
df$sample <- str_replace(df$sample, "Lok-75-Sample-4a-Ext-A12-Lib12A-Index1", "74a_B1_83_L4_Lok_75_Sample_4a_Ext_A12")
df$sample <- str_replace(df$sample, "Lok-75-Sample-4b-Ext-A25-Lib25A-Index1", "74a_B1_83_L4_Lok_75_Sample_4b_Ext_A25")
df$sample <- str_replace(df$sample, "Lok-75-Sample-4b-Ext-A26-Lib26A-Index1", "74a_B1_83_L4_Lok_75_Sample_4b_Ext_A26")
df$sample <- str_replace(df$sample, "12-1-47-Ext-26-Lib-26-Index1_C", "119_B3_116_L3_KapK_12_1_47_C")
df$sample <- str_replace(df$sample, "12-1-48-Ext-27-Lib-27-Index1_C", "119_B3_116_L3_KapK_12_1_48_C")
df$sample <- str_replace(df$sample, "12-1-49-Ext-28-Lib-28-Index1_C", "119_B3_116_L3_KapK_12_1_49_C")
df$sample <- str_replace(df$sample, "12-1-51-Ext-33-Lib-33-Index2_C", "119_B3_116_L4_KapK_12_1_51_C")
df$sample <- str_replace(df$sample, "12-1-52-Ext-34-Lib-34-Index1_C", "119_B3_116_L4_KapK_12_1_52_C")
df$sample <- str_replace(df$sample, "198A-Ext-55-Lib-55-Index1_C", "50_B3_127_L0_KapK_198A_C")
df$sample <- str_replace(df$sample, "199B-Ext-48-Lib-48-Index1_C", "50_B3_127_L0_KapK_199A_C")
df$sample <- str_replace(df$sample, "KapK-12-1-41-Ext-20-Lib-20-Index1_C", "119_B3_116_L1_KapK_12_1_41_C")
df$sample <- str_replace(df$sample, "KapK-12-1-42-Ext-21-Lib-21-Index1_C", "119_B3_116_L1_KapK_12_1_42_C")
df$sample <- str_replace(df$sample, "KapK-12-1-43-Ext-22-Lib-22-Index1_C", "119_B3_116_L1_KapK_12_1_43_C")
df$sample <- str_replace(df$sample, "KapK-12-1-45-Ext-24-Lib-24-Index1_C", "119_B3_116_L2_KapK_12_1_45_C")
df$sample <- str_replace(df$sample, "KapK-12-1-46-Ext-25-Lib-25-Index1_C", "119_B3_116_L2_KapK_12_1_46_C")
df$sample <- str_replace(df$sample, "200A-Ext-49-Lib-49-Index1_C", "50_B3_127_L0_KapK_200A_C")
df$sample <- str_replace(df$sample, "202B-Ext-38-Lib-38-Index2_C", "74b_B1_83_L1_KapK_202B_C")
df$sample <- str_replace(df$sample, "203A-Ext-42-Lib-42-Index1_C", "74b_B1_83_L1_KapK_203A_C")
df$sample <- str_replace(df$sample, "203B-Ext-43-Lib-43-Index1_C", "74b_B1_83_L1_KapK_203B_C")
df$sample <- str_replace(df$sample, "205A-Ext-51-Lib-51-Index1_C", "74b_B1_83_L2_KapK_205A_C")
df$sample <- str_replace(df$sample, "205B-Ext-52-Lib-52-Index2_C", "74b_B1_83_L2_KapK_205B_C")
df$sample <- str_replace(df$sample, "205C-Ext-53-Lib-53-Index1_C", "74b_B1_83_L2_KapK_205C_C")
df$sample <- str_replace(df$sample, "205D-Ext-54-Lib-54-Index2_C", "74b_B1_83_L2_KapK_205D_C")
df$sample <- str_replace(df$sample, "KapK-12-1-24-Ext-1-Lib-1-Index2_C", "119_B3_116_L0_KapK_12_1_24_C")
df$sample <- str_replace(df$sample, "KapK-12-1-25-Ext-2-Lib-2-Index1_C", "119_B3_116_L0_KapK_12_1_25_C")
df$sample <- str_replace(df$sample, "KapK-12-1-27-Ext-4-Lib-4-Index2_C", "119_B3_116_L0_KapK_12_1_27_C")
df$sample <- str_replace(df$sample, "KapK-12-1-29-Ext-6-Lib-6-Index2_C", "119_B3_116_L0_KapK_12_1_29_C")
df$sample <- str_replace(df$sample, "Lok-75-Sample-1-Ext-58-Lib-58-Index1_C", "74a_B1_83_L1_Lok_75_Sample_1_Ext_58_C")
df$sample <- str_replace(df$sample, "Lok-75-Sample-3-Ext-59-Lib-59-Index2_C", "74a_B1_83_L3_Lok_75_Sample_3_Ext_59_C")
df$sample <- str_replace(df$sample, "KapK-12-1-31-Ext-8-Lib-8-Index2_C", "69_B2_97_L0_KapK_12_1_31_C")
df$sample <- str_replace(df$sample, "KapK-12-1-33-Ext-10-Lib-10-Index2_C", "69_B2_97_L0_KapK_12_1_33_C")
df$sample <- str_replace(df$sample, "KapK-12-1-34-Ext-11-Lib-11-Index2_C", "69_B2_100_L0_KapK_12_1_34_C")
df$sample <- str_replace(df$sample, "KapK-12-1-35-Ext-12-Lib-12-Index2_C", "69_B2_100_L0_KapK_12_1_35_C")
df$sample <- str_replace(df$sample, "KapK-12-1-36-Ext-13-Lib-13-Index2_C", "69_B2_100_L0_KapK_12_1_36_C")
df$sample <- str_replace(df$sample, "KapK-12-1-37-Ext-17-Lib-17-Index1_C", "69_B2_103_L0_KapK_12_1_37_C")
df$sample <- str_replace(df$sample, "KapK-12-1-38-Ext-18-Lib-18-Index1_C", "69_B2_103_L0_KapK_12_1_38_C")
df$sample <- str_replace(df$sample, "KapK-12-1-39-Ext-19-Lib-19-Index1_C", "69_B2_103_L0_KapK_12_1_39_C")
unique(df$sample)

# Minimum amount of damage (filtered for)
DamMin = 0.250
#Lambda Likelihood Ratio
LR = 1.5
# Minimum reads for parsing taxa
MinRead = 5
# Minimum mean read length
MinLength = 30
###### subsetting the table using grepl and filter, parameters you need to set and possible add more
df1 <- df %>% filter(D_max > DamMin, N_reads >= MinRead, mean_L > MinLength, lambda_LR  > LR,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))
unique(df1$sample)
unique(df$sample)

# dataframe for all plant taxa with no damage filtering, but with 500 reads or more, for DNA damage histogram
DamMin2 = 0.00
LR2 = 0
MinRead2 = 500
df2 <- df %>% filter(D_max > DamMin2, N_reads >= MinRead2, mean_L > MinLength, lambda_LR  > LR2,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))

mean(df2$D_max)
median(df2$D_max)

#saved as 4x6 figure
histogram(df2$D_max, col="darkmagenta", breaks=200, xlab="D_max")

###### converting the long table to wide format
# wide table with N_reads
df1a <- as.data.frame(df1)
data_wide <- dcast(df1a, tax_name ~ sample, value.var="N_reads")

# optional you can make the table with any value here wide table with Dmax
#data_wide2 <- dcast(subset5, tax_name ~ sample, value.var="D_max")

n <- ncol(data_wide)
dw1 <- data_wide[,2:n]
rownames(dw1) <- data_wide$tax_name
dw1[is.na(dw1)] <- 0 #set all NAs as zeros

######### setting thresholds for minimum nuber of reads per taxa and per sample to be parsed, plots will show you later how many you remove they can be adjusted accordingly.  
colSums(dw1) # this prints sum of each sample (or column) for the subsetted wide table
col <- median.default(colSums(dw1))/2 # threshold1 removes all samples with less than half the median of the sum os reads in all samples
row <- median.default(rowSums(dw1))/2 # threshold2 removes all taxa with less than half the median of the sum os reads to all taxa 

drop <- colSums(dw1) > col # applies the threshold1
dd1 <- dw1[c(drop)] # applies the threshold1
dd2 <- as.data.frame(t(dd1)) # tranposes the dataframe  
drop2 <- colSums(dd2) > row # applies the threshold2
dd3 <- dd2[c(drop2)] # applies the threshold2
dd4 <- as.data.frame(t(dd3)) # transposes the dataframe (df) back to original

#hist(as.matrix(dw1), breaks = 100000) # plots histogram of the original df with taxa and their number of reads
#hist(as.matrix(dd4), breaks = 100000) # plots histogram of the filtered taxa and their number of reads

# adding new column with sum of reads
b1 <- cbind(dd4, NReads = rowSums(dd4))
i=ncol(b1)-1
b2 <- cbind(b1, Nreplicated = rowSums(b1[,1:i] > 0))

# how many times must a taxa appear in different samples to be considered 'replicated', considers all samples as different. 
b2 <- b2[b2$Nreplicated >= 3,] # applies the replication threshold

#prints sum of sampels and taxa and samples names
colSums(b2)
rowSums(b2)
colnames(b2)

#test, if this one fails it might have text in the number of reads coloumn
b2[is.na(b2)]=0 


##################### create percentage table #######################
# calculated the number of columns and subtracts the sum of all reads coloumn, remove the -1 and you will plot the Nreads column on later plots
i=ncol(b2)-1
b3=as.matrix(b2[,seq(1,i)])  # change number of coloumns to the total of mydata

b4 <- prop.table(data.matrix(b3), margin=2)*100 # makes proportion table, needs 2 margins e.g. header and 1st row names
colSums(prop.table(b4, margin=2)*100) # prints sum of column, should give 100 for each 

b4[is.nan(b4)]=0 
############### setting a 3rd threshold for how big a taxa needs of percentage to be parsed
j <- median.default(rowSums(b4))*3

# prints results from logical statement, TRUE are parsed for plotting FALSE are removed
rowSums(b4) > j
b5 <- b4[apply(b4[,1:i], MARGIN = 1, function(x) any(x > j)), ] # applies the threshold

############## removing only rows if no values are above x in all row boxes OBS alternative threshold
# tmp2 <- tmp[apply(b3[,1:i], MARGIN = 1, function(x) all(x > 0.01)), ]


### define the colors within 4 zones for the Hclust heatmaps
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
breaks = seq(0, max(b5), length.out=100)
gradient1 = colorpanel( sum( breaks[-1]<=45 ), "white", "green" )
gradient2 = colorpanel( sum( breaks[-1]>45 & breaks[-1]<=65 ), "green", "blue" )
gradient3 = colorpanel( sum( breaks[-1]>65 & breaks[-1]<=95 ), "blue", "darkblue" )
gradient4 = colorpanel( sum( breaks[-1]>95 ), "darkblue", "black" )

hm.colors = c(gradient1, gradient2, gradient3, gradient4)

# can output the dataframe to comma seperatefile if needed
now<-format(Sys.time(), "%d-%m-%Y")
wd <- getwd()
csvFileName <- paste("/",wd,now,"KapKplantGenusMetaDMG6out.csv",sep="")
write.csv(b5, file=csvFileName)

## removing Nreads column
z <- ncol(b5)
y <- z-1

# makes b5 to long table for ggplot heatmap below
b6 <- melt(b5[,1:y])
colnames(b6) <- c("Taxa", "Sample", "Percentage")

## nMDS 
nmds <- metaMDS(b5[,1:y],trymax = 200)

### converting nmds to long format for ggplot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- colnames(nmds)  #  add the grp variable created earlier
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores, 50)  #look at the data

p0 <- ggplot() + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0) +  # add the site labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, fill="grey"),size=0.5,alpha=0.5) +  # add the site labels
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,colour=species),size=3) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species, colour=species), vjust=-1) +  # add the species labels
  ggtitle("NMDS ordination of filtered samples (all plots following have filtered data)") +
  labs(col = "Samples", fill="Driving taxa") +
  theme_bw()

p0 # plots nmds plot in your GUI

## filter only taxa and genus used and plot a damageplot prior to filtering 
taxa1 <- rownames(dd4)
dm1 <- df1[df1$tax_name %in% taxa1,]

## filter only taxa and genus used and plotted in heatmap (after filtering)
taxa2 <- unique(b6$Taxa)
dm2 <- df1[df1$tax_name %in% taxa2,]

# making table with samples parsed only
sample2 <- unique(b6$Sample)
dfx2 <- dm2[dm2$sample %in% sample2,]

unique(dfx2$sample)
unique(dfx2$tax_name)
rawData <- dcast(dfx2, sample ~ tax_name, value.var="N_reads")
rawData2 <- t(as.data.frame(rawData))
rawData2[is.na(rawData)]=0 
write.csv(rawData2, file="KapKmetaDMG6PlantGenusRawReadCountOut.csv")

mean(dfx2$D_max)

## making table with taxa parsed after filters, but with raw reads. 
data_wideFiltered <- data_wide[data_wide$tax_name %in% taxa2,]
data_wideFiltered[is.na(data_wideFiltered)]=0

ncol(data_wideFiltered)
drop4 <- colSums(data_wideFiltered[,2:ncol(data_wideFiltered)]) > 0
drop4 <- colSums(data_wideFiltered[,2:ncol(data_wideFiltered)]) > 0
tmp1 <- data_wideFiltered[,2:ncol(data_wideFiltered)]
data_wideFiltered2 <- tmp1[c(drop4)] # applies the threshold2
rownames(data_wideFiltered2) <- data_wideFiltered$tax_name

write.csv(data_wideFiltered2, file="KapKmetaDMG6PlantGenusReadCountOut.csv")


pdf("KapKPlantGenusMetaDMG6out.pdf", height = 10,width = 18)
text1=paste("This pdf contains a series of plots of the Plant Genus found by DNA from KapK 
Below are the general thresholds set, which have been set to minimize random mapping noise.

Threshold values:", " 
Minimum DNA damage (Dmax) = ",DamMin,"
Minimum likelihood ratio (LR) = ",LR,"
Minimum mean read length = ",MinLength,"
Threshold 1 minimum reads per taxa = ",row,", 
Threshold 2 minimum reads per sample = ", col,", 
Threshold 3 minimum percentage per taxa = ", j, sep = "")

### test plotting
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,text1, pos=4)


df77 = data.frame(Samples=colnames(dw1),Total_reads=colSums(dw1))
ggplot(data=df77, aes(x=Samples, y=Total_reads)) +
  geom_bar(stat="identity")  +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))


# Plotting length per taxa as boxplots
ggplot(dm2, aes(x=mean_L, y=tax_name)) + geom_boxplot() + theme_test() + ggtitle("Mean read lengths for each sample per genus") +
  ylab("Taxa") + xlab("Mean Length (Bp)") 
# Plotting Dmax per taxa as boxplots
ggplot(dm2, aes(x=D_max, y=tax_name)) + geom_boxplot() + theme_test() + ggtitle("Average DNA damage on 1st position of the read (Dmax) per taxa") +
  ylab("Taxa") + xlab("Dmax")  + xlim(0.15,0.85) 
# Plotting Dmax per sample as boxplots
ggplot(dm2, aes(x=D_max, y=sample)) + geom_boxplot() + theme_test() + ggtitle("Average DNA damage on 1st position of the read (Dmax) per sample") +
  ylab("Taxa") + xlab("Dmax")  + xlim(0.15,0.85) 
# Plotting length per sample as boxplots
ggplot(dm2, aes(x=mean_L, y=sample)) + geom_boxplot() + theme_test() + ggtitle("Mean read lengths for each genus per sample") +
  ylab("Taxa") + xlab("Mean Length (Bp)") 

## uncategorised plant heatmap, not part of the publication but are further down sorted by growth form and plotted.
p2 <- ggplot(b6, aes(x=Taxa, y=Sample, fill=log10(Percentage))) +   geom_tile(colour="lightgrey") + 
  theme_minimal() + scale_fill_gradient(low="white", high="darkgreen") + scale_y_discrete(limits=rev)

p2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Logbase 10 transformed percentage of taxa plotted as heatmap (highlights less abundant taxa)") +
  xlab("Taxa") + ylab("Samples") + labs(fill = "Logbase10 transformed \npercentage")

dev.off()


### final plots for publication ##############################################################################3

##### merging growth classes with genus tables, growth classes is a list of genus, growth form, the rank you want each growthform to come, presence absence in the fossil record.
#plantGenusDNA_SF050921 <- read_excel("~/Dropbox/Kap_København/Macrofossil data/DNA_Fossil_comparison/oldversions/plantGenusDNA_SF050921.xlsx", sheet = "growth_merge")
#sb <- plantGenusDNA_SF050921[,2:5]

colnames(b6) <- c("Taxa", "Sample", "Percentage")
b6dw <- dcast(b6, Taxa ~ Sample, value.var="Percentage")
growth_IGA <- read_excel("/data/growth_form_table.xlsx", sheet = "GrowthFormTable5fe2022_reduced")

KapK_df2 <- merge.data.frame(b6dw , growth_IGA, by=c("Taxa"),all=TRUE)

## saving the plant table sorted by growth form
newdata2 <- KapK_df2[t(order(KapK_df2$Growth_form_order, decreasing = TRUE)),]
write.csv(newdata2, file = 'KapKmetaDMG6_PlantGenusDNAgrowthForm_added.csv')

#cutting away growth information, transposing data and converting to long table for plotting
v <- ncol(newdata2)-5
dd1 <- newdata2[,2:v]
rownames(dd1) <- newdata2$Taxa
sapply(dd1, class) 
colSums(dd1)

dd5 <- melt(as.matrix(dd1))

pdf("KapKmetaDMG6_MainPlantFigLog10_transposed.pdf", height = 16.5,width = 9)
p <- ggplot(dd5, aes(y=Var1, x=Var2, fill=log10(value))) +   geom_tile(colour="lightgrey") + 
  theme_minimal() + scale_fill_gradient(low="white", high="darkgreen") 

p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust =1))
dev.off()


##### final DNA damage plots
pdf("KapKmetaDMG6_PlantGenusDNA14x7DmaxPerTaxa.pdf", height = 14,width = 7)
# Plotting Dmax per taxa as boxplots
ggplot(dm2, aes(x=D_max, y=tax_name)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Average DNA damage on 1st position of the read (Dmax) per taxa") +
  ylab("Taxa") + xlab("Dmax")  + xlim(0.15,0.85) 
dev.off()
pdf("KapKmetaDMG6_PlantGenusDNA14x7DmaxPerSample.pdf", height = 14,width = 7)
# Plotting Dmax per sample as boxplots
ggplot(dm2, aes(x=D_max, y=sample)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Average DNA damage on 1st position of the read (Dmax) per sample") +
  ylab("Taxa") + xlab("Dmax")  + xlim(0.15,0.85) 
dev.off()

##### final Read lenght plots
pdf("KapKmetaDMG6_PlantGenusDNA14x7ReadlengthPerSample.pdf", height = 14,width = 7)
ggplot(dm2, aes(x=mean_L, y=sample)) + geom_boxplot() + theme() + ggtitle("Mean read lengths for each genus per sample") +
  ylab("Taxa") + xlab("Mean read length (bp)") 
dev.off()

pdf("KapKmetaDMG6_PlantGenusDNA14x7ReadlengthPerTaxa.pdf", height = 14,width = 7)
# Plotting length per taxa as boxplots
ggplot(dm2, aes(x=mean_L, y=tax_name)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Boxplot of DNA length in basepairs (mean-L) for each sample (filtered)") +
  xlab("Mean read length (bp)") + ylab("Mean Length (Bp)")
dev.off()
