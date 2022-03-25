library(tidyverse)
library(janitor)
library(adephylo)
library(ape)
library(phytools)
library(ggpubr)
library(readxl)
library(ggthemr)

source("libs/lib.R")
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

site_order <- tibble(
  site = c("119", "50", "69", "75", "74", "nosite"),
  site_rnk = 1:6
)


# Read Delmont et al. SMAGS data ------------------------------------------

smags_metadata <- read_tsv("data/map-smags/SMAGS-anvio/metadata.txt") %>%
  select(Genome_Id, contains("Best_taxonomy"), contains("119_Stations"), Taxa_Super_Groups) %>%
  mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))

# Mapping stats
# We will keep those SMAGS with:
# more than 500 reads mapped
# the ratio between the observed and expected breadth ratio >= 0.75
# the coverage evenness >= 0.75
# and a mean read ANI >= 93%

nreads <- 500
oe_bratio <- 0.75
cov_e <- 0.75
ani <- 93
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
mtdmg_w1_all <- read_tsv(file = "data/map-smags/smags-mdmg.summary.1.tsv.gz") %>%
  separate(tax_name, into = c("name", "Genome_Id"), sep = " ", remove = FALSE)

mtdmg_w1_filt_all <- mtdmg_w1_all %>%
  filter(N_reads > nreads) %>%
  filter(Bayesian_z >= 1.5, Bayesian_D_max >= 0.25) %>%
  inner_join(kapk_cdata %>% select(label, member_unit)) %>%
  arrange(desc(Bayesian_D_max)) %>%
  mutate(Genome_Id = gsub("Metagenome_centric_SAG_", "", Genome_Id))


mtdmg_w1 <- read_tsv(file = "data/map-smags/smags-mdmg.summary.1.tsv.gz") %>%
  filter(tax_rank == "species") %>%
  separate(tax_name, into = c("name", "Genome_Id"), sep = " ", remove = FALSE)

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

# We kept a total of 24 SMAGS with both filtering steps
smags_stats %>%
  pull(Genome_Id) %>%
  unique() %>%
  length()


# Some supp figure with number of reads and ANI ---------------------------

spg_colors <- c("#44AD75", "#4EB1DC", "#E86DB7")
names(spg_colors) <- c("Opisthokonta", "Stramenopiles", "Archaeplastida")
smags_stats %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt) %>%
  ggplot(aes(read_ani_mean, n_reads, fill = Taxa_Super_Groups, size = Bayesian_D_max)) +
  geom_point(shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 7), limits = c(0.25, 0.6)) +
  xlab("Mean read ANI") +
  ylab("Number of reads") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_fill_manual(values = spg_colors)



# Preparing the data for anvi'o -------------------------------------------

# We read the phylogenomic tree
smags_tree <- ape::read.tree("data/map-smags/SMAGS-anvio/Conca_allMETDB_allSMAGs_RNAP_200aa_no_duplicate_NUCLEAR_noRhodophyta.fa.treefile")

# Read taxonomic annotations for the SMAGS
smags_tax <- read_tsv("data/map-smags/SMAGS-anvio/TARA_SMAGs.txt") %>%
  select(accession, taxon) %>%
  tidyr::separate(col = taxon, sep = ";", into = c("domain", "phylum", "class", "order", "family", "genus", "species")) %>%
  rename(Genome_Id = accession)

# Let;s do a quick fix of the taxonomies
smags_tax %>%
  left_join(smags_metadata %>%
    select(Genome_Id, Taxa_Super_Groups)) %>%
  mutate(genus = ifelse(genus == "g__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(genus, "_", Taxa_Super_Groups), genus)) %>%
  mutate(family = ifelse(family == "f__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(family, "_", Taxa_Super_Groups), family)) %>%
  mutate(order = ifelse(order == "o__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(order, "_", Taxa_Super_Groups), order)) %>%
  mutate(class = ifelse(class == "c__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(class, "_", Taxa_Super_Groups), class)) %>%
  unite(domain, phylum, class, order, family, genus, species, col = "tax_path", sep = ";") %>%
  select(-Taxa_Super_Groups) %>%
  write_tsv(file = "data/map-smags/SMAGS-anvio/smags-tax-fixed.tsv", col_names = F)

# Not all SMAGS are present in the tree because they didn't have RNApol, let's do a bit of a hack
# We will use a anchor the taxonomic information from Delmont's best taxonomy. We will first
# try to get as many in the lower ranks as we can. Starting from Genus to Kingdom
# For example:
# we will identify those SMAGS not present in the tree, get their genera and find out where in the tree they are
# then create a subtree using those tips and get the node that it closer to the root.
# We will place in this node the SMAGS that are missing

# Which genera have a taxonomy
miss_gen <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(!(Genome_Id %in% smags_tree$tip.label)) %>%
  select(Best_taxonomy_GENRE, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_GENRE)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_GENRE)

# Which families have a taxonomy
miss_fam <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id)
  ) %>%
  select(Best_taxonomy_FAMILY, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_FAMILY)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_FAMILY)


# Which orders have a taxonomy
miss_order <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id)
  ) %>%
  select(Best_taxonomy_ORDER, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_ORDER)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_ORDER)


# Which classes have a taxonomy
miss_class <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id),
    !(Genome_Id %in% miss_order$Genome_Id)
  ) %>%
  select(Best_taxonomy_CLASS, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_CLASS)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_CLASS)

# Which phyla have a taxonomy
miss_phyl <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id),
    !(Genome_Id %in% miss_order$Genome_Id),
    !(Genome_Id %in% miss_class$Genome_Id)
  ) %>%
  select(Best_taxonomy_PHYLUM, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_PHYLUM)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_PHYLUM)

# Which Kingdoms have a taxonomy
miss_king <- smags_stats %>%
  select(Genome_Id, label) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id),
    !(Genome_Id %in% miss_order$Genome_Id),
    !(Genome_Id %in% miss_class$Genome_Id),
    !(Genome_Id %in% miss_phyl$Genome_Id)
  ) %>%
  select(Best_taxonomy_KINGDON, Genome_Id) %>%
  filter(!is.na(Best_taxonomy_KINGDON)) %>%
  distinct() %>%
  rename(tax = Best_taxonomy_KINGDON)

miss_smags <- bind_rows(miss_gen, miss_fam, miss_order, miss_class, miss_phyl, miss_king)

# Which ones are in the tree
in_tree <- smags_stats %>%
  select(Genome_Id) %>%
  distinct() %>%
  inner_join(smags_metadata) %>%
  inner_join(mtdmg_w1_filt %>% select(Genome_Id, label)) %>%
  distinct() %>%
  filter((Genome_Id %in% smags_tree$tip.label))

# Find which ones to drop
to_drop <- setdiff(smags_tree$tip.label, unique(smags_stats$Genome_Id))

# We drop as many
length(to_drop)

# Filter tree
tree_filt <- drop.tip(smags_tree, tip = to_drop)
# tree_filt <- smags_tree

# Get tip's info
get_tip_info <- function(X, df, phy_tree) {
  cat(paste(X, "\n"))
  tips <- df %>%
    filter(tax == X)
  to_drop <- setdiff(phy_tree$tip.label, unique(tips$Genome_Id))
  phy_tree_filt <- drop.tip(phy_tree, tip = to_drop)
  if (!is.null(phy_tree_filt)) {
    if (phy_tree_filt$Nnode == 1) {
      node_n <- which(phy_tree$tip.label == phy_tree_filt$tip)
      tibble(tip = phy_tree_filt$tip.label, dist = phy_tree_filt$edge.length, tax = X) %>% mutate(node = node_n)
    } else {
      nodes <- distRoot(phy_tree_filt, tips = "all") %>%
        sort(decreasing = F) %>%
        head(1) %>%
        enframe(name = "tip", value = "dist") %>%
        mutate(tax = X)
      node_n <- which(phy_tree$tip.label == nodes$tip)
      nodes %>% mutate(node = node_n)
    }
  }
}

# get tips based on the ranks we defined before
df <- smags_metadata %>%
  select(Genome_Id, Best_taxonomy_GENRE) %>%
  rename(tax = Best_taxonomy_GENRE) %>%
  distinct()

gen_tips <- map_dfr(miss_gen %>% select(tax) %>% distinct() %>% pull(tax),
  get_tip_info,
  df = df, phy_tree = smags_tree
) %>%
  inner_join(miss_gen)

# Combine
df_all <- bind_rows(gen_tips)

# Create a new tree where we plug the missing genomes closer to the existing ranks
new_tree <- smags_tree
for (i in 1:nrow(df_all)) {
  new_tree <- bind.tip(new_tree,
    tip.label = df_all[i, ]$Genome_Id,
    where = df_all[i, ]$node,
    edge.length = df_all[i, ]$dist
  )
}


to_drop <- setdiff(new_tree$tip.label, unique(smags_stats$Genome_Id))
smags_tree_filt <- drop.tip(new_tree, tip = to_drop)

# Write final tree
write.tree(smags_tree_filt, file = "data/map-smags/SMAGS-anvio/smags-tree.tree")

# Keep only those with a sigma and model
mtdmg_w1_filt_anv <- mtdmg_w1_filt %>%
  dplyr::select(Genome_Id, Bayesian_D_max, Bayesian_z, label)

# Prepare data to be plotted in anvi'o
# Get damage values
dmg <- smags_stats %>%
  inner_join(mtdmg_w1_filt) %>%
  filter(Genome_Id %in% smags_tree_filt$tip.label) %>%
  mutate(Bayesian_D_max = ifelse(is.na(Bayesian_D_max), 0, Bayesian_D_max)) %>%
  mutate(member_unit = paste0(member_unit, "_Bayesian_D_max")) %>%
  select(Genome_Id, member_unit, Bayesian_D_max) %>%
  group_by(Genome_Id, member_unit) %>%
  summarise(avg_D_max = mean(Bayesian_D_max)) %>%
  pivot_wider(names_from = member_unit, values_from = avg_D_max, values_fill = 0)

n_reads <- smags_stats %>%
  inner_join(mtdmg_w1_filt) %>%
  select(Genome_Id, member_unit, n_reads) %>%
  group_by(Genome_Id, member_unit) %>%
  summarise(avg_mapped_reads = sum(n_reads)) %>%
  mutate(member_unit = paste0(member_unit, "_nreads")) %>%
  filter(Genome_Id %in% smags_tree_filt$tip.label) %>%
  pivot_wider(names_from = member_unit, values_from = avg_mapped_reads, values_fill = 0)

exp_br <- smags_stats %>%
  inner_join(mtdmg_w1_filt) %>%
  select(Genome_Id, member_unit, breadth_exp_ratio) %>%
  group_by(Genome_Id, member_unit) %>%
  summarise(avg_breadth_exp_ratio = mean(breadth_exp_ratio)) %>%
  mutate(member_unit = paste0(member_unit, "_breadth_exp_ratio")) %>%
  filter(Genome_Id %in% smags_tree_filt$tip.label) %>%
  pivot_wider(names_from = member_unit, values_from = avg_breadth_exp_ratio, values_fill = 0)

cosm_reg <- c(
  "Arctic_119_Stations", "Atlantic_119_Stations", "Indian_119_Stations", "Mediterranean_119_Stations",
  "Red Sea_119_Stations", "Pacific_119_Stations", "Southern_119_Stations"
)

regions <- smags_metadata %>%
  filter(Genome_Id %in% smags_stats$Genome_Id) %>%
  # filter(Genome_Id %in% smags_order) %>%
  select(Genome_Id, all_of(cosm_reg)) %>%
  janitor::clean_names() %>%
  rename(Genome_Id = genome_id)

names(regions) <- c("Genome_Id", paste0("bars!", names(regions)[2:length(names(regions))]))

n_reads %>%
  # inner_join(exp_br) %>%
  inner_join(dmg) %>%
  inner_join(smags_metadata %>% select(Genome_Id, Taxa_Super_Groups)) %>%
  inner_join(regions) %>%
  filter(Genome_Id %in% smags_tree_filt$tip.label) %>%
  write_tsv("data/map-smags/SMAGS-anvio/smags-cov-dmg-ebr.tsv")


# DECAY fit ---------------------------------------------------------------
# Let's plot the damage profile for each SMAG in each sample
# Get the taxa supergroups
tax_superg <- smags_stats %>%
  inner_join(smags_metadata) %>%
  select(Genome_Id, Taxa_Super_Groups) %>%
  mutate(species = Genome_Id)


tax_g_list <- smags_stats %>%
  inner_join(smags_metadata) %>%
  pull(Taxa_Super_Groups) %>%
  unique()


mtdmg_w1_stats <- mtdmg_w1_filt %>%
  inner_join(smags_stats)


ggthemr(layout = "scientific", palette = "fresh")
purrr::map(tax_g_list, function(X, nrank) {
  sel_tax <- mtdmg_w1_stats %>%
    filter(Genome_Id %in% (tax_superg %>% filter(Taxa_Super_Groups == X) %>% pull(Genome_Id))) %>%
    filter(tax_rank == nrank) %>%
    select(tax_name, label) %>%
    distinct() %>%
    arrange(tax_name)
  if (nrow(sel_tax) > 0) {
    n_readsa <- mtdmg_w1_stats %>%
      inner_join(sel_tax) %>%
      filter(tax_rank == nrank) %>%
      pull(N_reads) %>%
      sum()
    ggpubr::ggarrange(plotlist = list(
      get_dmg_decay_fit(df = mtdmg_w1_stats %>% inner_join(sel_tax) %>% filter(tax_rank == nrank), orient = "fwd") +
        ggtitle(paste0(X, " ", n_readsa, " FWD")),
      get_dmg_decay_fit(df = mtdmg_w1_stats %>% inner_join(sel_tax) %>% filter(tax_rank == nrank), orient = "rev") +
        ggtitle(paste0(X, " ", n_readsa, " REV"))
    ), align = "hv")
    ggsave(paste0("figures/", X, "-dmg.pdf"), plot = last_plot(), width = 8, height = 4)
  }
}, nrank = "species")

# DECAY fit end -----------------------------------------------------------
