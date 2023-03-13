library(tidyverse)
library(janitor)
library(R.utils)

# Download and load files
options(timeout=300)
download.file("https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz", "data/bac120_metadata_r207.tar.gz")
download.file("https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz", "data/ar53_metadata_r207.tar.gz")

# Transform compression
untar("data/bac120_metadata_r207.tar.gz")
untar("data/ar53_metadata_r207.tar.gz")

R.utils::gzip("data/bac120_metadata_r207.tsv")
R.utils::gzip("data/ar53_metadata_r207.tsv")

# Load tables
bac <- read_delim("data/bac120_metadata_r207.tsv.gz")
arch <- read_delim("data/ar53_metadata_r207.tsv.gz")

gtdb207 <- rbind(bac, arch)

# Get/plot genome stats (>99 complete, <5 contam, >50 contigs)
gtdb207_filt <- gtdb207 %>%
  filter(checkm_completeness > 99 & checkm_contamination < 5 & scaffold_count > 50) %>%
  separate(., col = gtdb_taxonomy, sep = ";", into = c("domain", "phylum", "class", "order", "family", "genus", "species"))

z<-gtdb207_filt %>%
  group_by(species) %>%
  tally()

#Choose only 1 representative per species
gtdb207_filt_drep_species <- gtdb207_filt %>%
  group_by(species) %>%
  sample_n(1)

# Plot
gtdb207_filt_drep_species %>%
  group_by(phylum) %>%
  tally(sort = T) %>%
  ggplot(aes(x = fct_reorder(phylum, n), y = n)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("figures/selection_phyla_representation.pdf", width = 7, heigh = 7)
