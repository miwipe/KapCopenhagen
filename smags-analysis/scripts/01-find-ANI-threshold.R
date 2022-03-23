library(tidyverse)
library(PCAtools)
# Some supp figure with number of reads and ANI ---------------------------

# Read KapK metadata ------------------------------------------------------
kapk_cdata <- readxl::read_xlsx("data/cdata/KapK_samples-20210702.xlsx") %>%
  clean_names() %>%
  select(member_unit, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, figure_names, site)

kapk_cdata <- kapk_cdata %>%
  mutate(
    label = basename(collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads),
    label = gsub("\\..*", "", label)
  ) %>%
  distinct() %>%
  select(label, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, member_unit, figure_names, site) %>%
  rename(label_orig = label, label = collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum) %>%
  mutate(
    site = as.character(site),
    site = ifelse(site == "NA", "nosite", site),
    figure_names = paste0(figure_names, "-", substr(label, 1, 5))
  )

# Explore the number read decay as function of the mean read ANI values
anis <- seq(90, 100, 1)
mtdmg_w1 <- read_tsv(file = "data/map-smags/smags-mdmg.summary.1.tsv.gz") %>%
  filter(tax_rank == "species") %>%
  separate(tax_name, into = c("name", "Genome_Id"), sep = " ", remove = FALSE)

# Set some filtering values
nreads <- 500
oe_bratio <- 0.75
cov_e <- 0.75
anis_n <- map_dfr(anis, function(ani) {
  smags_stats <- read_tsv(file = "data/map-smags/smags-mapping.summary.tsv.gz") %>%
    filter(breadth_exp_ratio >= oe_bratio, n_reads >= nreads, cov_evenness >= cov_e, read_ani_mean >= ani) %>%
    select(reference, n_reads, label, breadth_exp_ratio, cov_evenness, read_ani_mean, edit_distances) %>%
    rename(Genome_Id = reference) %>%
    inner_join(kapk_cdata) %>%
    mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))
  # Get metaDMG estimates ---------------------------------------------------
  # Here as each SMAG is equivalent to one species we used the SMAG name as
  # species name to get genome-wise damage estimates
  # Read damage estimates
  # we use weight type 1: mismatch matrix is weighted or rescaled according to number of alignmetns
  # Read metaDMG estimates and filter
  mtdmg_w1_filt <- mtdmg_w1 %>%
    filter(N_reads > nreads) %>%
    filter(Bayesian_z >= 1.5, Bayesian_D_max >= 0.25) %>%
    inner_join(kapk_cdata %>% select(label, member_unit)) %>%
    arrange(desc(Bayesian_D_max)) %>%
    mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))
  # Which smags are damaged
  dmg_smags <- mtdmg_w1_filt %>%
    select(Genome_Id) %>%
    distinct() %>%
    pull(Genome_Id)
  # Get smags that have a good sigma
  smags_stats <- smags_stats %>%
    inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
    distinct()
  # We kept a total of 317 SMAGS with both filtering steps
  tibble(
    ani = ani,
    n_mags = smags_stats %>%
      pull(Genome_Id) %>%
      unique() %>%
      length(),
    n_reads = smags_stats %>% pull(n_reads) %>% sum()
  )
})

# identify the elbow point
elbow_point <- anis_n %>% filter(n_reads == anis_n$n_reads[findElbowPoint(anis_n$n_reads)])

ggthemr::ggthemr(layout = "scientific", palette = "fresh")

anis_n %>%
  ggplot(aes(ani, n_reads), group = 1) +
  geom_line() +
  geom_point(size = 3, shape = 21, color = "black", fill = "red") +
  gghighlight::gghighlight(ani == elbow_point$ani, unhighlighted_params = list(fill = "grey")) +
  scale_x_continuous(breaks = anis) +
  scale_y_continuous(labels = scales::comma) +
  xlab("Mean read ANI (%)") +
  ylab("Number of reads")


# Let's see the number of reads per each super group and ANI value

# Read the SMAGs metadata
smags_metadata <- read_tsv("data/map-smags/SMAGS-anvio/metadata.txt") %>%
  select(Genome_Id, contains("Best_taxonomy"), contains("119_Stations"), Taxa_Super_Groups) %>%
  mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))

# Filter metaDMG estimates
mtdmg_w1_filt <- mtdmg_w1 %>%
  filter(N_reads > nreads) %>%
  filter(Bayesian_z >= 1.5, Bayesian_D_max >= 0.25) %>%
  inner_join(kapk_cdata %>% select(label, member_unit)) %>%
  arrange(desc(Bayesian_D_max)) %>%
  mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))


# Get mapping stats from filterBAM
spg_colors <- c("#44AD75", "#4EB1DC", "#E86DB7")
names(spg_colors) <- c("Opisthokonta", "Stramenopiles", "Archaeplastida")
read_tsv(file = "data/map-smags/smags-mapping.summary.tsv.gz") %>%
  filter(breadth_exp_ratio >= oe_bratio, n_reads >= nreads, cov_evenness >= cov_e, read_ani_mean >= elbow_point$ani) %>%
  select(reference, n_reads, label, breadth_exp_ratio, cov_evenness, read_ani_mean, edit_distances) %>%
  rename(Genome_Id = reference) %>%
  inner_join(kapk_cdata) %>%
  mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id)) %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt) %>%
  ggplot(aes(read_ani_mean, n_reads, fill = Taxa_Super_Groups, size = Bayesian_D_max)) +
  geom_point(shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 7), limits = c(0.25, 0.6)) +
  xlab("Mean read ANI (%)") +
  ylab("Number of reads") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_fill_manual(values = spg_colors)
