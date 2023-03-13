#setwd("/Users/anttonalberdi/github/2023_genome_function_2")

library(tidyverse)
library(broom)
library(ggplot2)
library(Rtsne)

#Load DRAM products
DRAM_products <- read_tsv("data/20K_merged_product.tsv.gz") %>%
  mutate_at(vars(genome), ~ str_extract(., "[^_]*_[^_]*")) #simplify to assembly accession


#Load metadata
genome_metadata <- read_delim("data/20K_genomes_metadata.tsv.gz") %>%
  rename(genome=accession) %>%
  mutate_at(vars(genome), ~ str_replace(., "GB_","")) %>%
  mutate_at(vars(genome), ~ str_replace(., "RS_",""))

#Merge
DRAM_products_metadata <- DRAM_products %>%
  left_join(genome_metadata, by="genome") %>%
  filter(!is.na(coding_density))  #filter genomes with metadata available

#####
# PCA
#####

#Run principal components analysis (PCA)
DRAM_products_pc <- DRAM_products_metadata %>%
  select(matches("M[0-9][0-9][0-9][0-9][0-9]")) %>% prcomp()

#Plot PCA
DRAM_products_pc %>%
    augment(DRAM_products_metadata) %>%
    ggplot(aes(.fittedPC1, .fittedPC2, color=phylum)) +
    geom_point(size = 1) +
    theme_minimal() +
    theme(legend.position = "none")

#####
# tSNE
#####

DRAM_products_tsne <- DRAM_products_metadata %>%
  select(matches("M[0-9][0-9][0-9][0-9][0-9]")) %>%
  Rtsne(X=., dims = 2, check_duplicates = FALSE)

DRAM_products_tsne$Y %>%
  as_tibble() %>%
  mutate(genome=DRAM_products_metadata$genome) %>%
  left_join(genome_metadata, by="genome") %>%
  rename(tSNE1="V1", tSNE2="V2") %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = phylum))+
    geom_point(size=0.5, shape=16, alpha=0.8) +
    theme_minimal() +
    theme(legend.position = "none")

#####
# Selection through clustering
#####
#https://stackoverflow.com/questions/73541633/how-can-i-subset-data-with-r-to-maximize-diversity-across-each-of-multiple-colum

###### First stage selection based
###### on functional features

k=5000 #number of desired clusters

#Calculate distances
DRAM_products_dist <- DRAM_products_metadata %>%
  select(matches("M[0-9][0-9][0-9][0-9][0-9]")) %>%
  dist()

#Divide into k clusters using hierarchical clustering
DRAM_products_clust <- DRAM_products_dist %>%
  hclust(., method = "complete") %>%
  cutree(hc, k = k)

# Select a single representative from each group
DRAM_products_selected <- DRAM_products_metadata %>%
    mutate(cluster = DRAM_products_clust) %>%
    group_by(cluster) %>%
    slice_sample(n = 1) %>%
    select(genome,cluster)

#Add whether selected to metadata
DRAM_products_metadata_cluster <- DRAM_products_metadata %>%
  left_join(DRAM_products_selected,by="genome") %>%
  relocate(cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), 0, 1))

###### Second stage selection based
###### on genomic features

k=1000

DRAM_products_metadata_cluster_dist <- DRAM_products_metadata_cluster %>%
  filter(cluster == 1) %>%
  select(coding_density,gc_percentage,genome_size,protein_count)  %>%
  mutate_all(~(scale(.) %>% as.vector)) %>%
  dist()

DRAM_products_metadata_cluster_clust <- DRAM_products_metadata_cluster_dist %>%
  hclust(., method = "complete") %>%
  cutree(hc, k = k)

# Select a single representative from each group
DRAM_products_metadata_cluster_selected <- DRAM_products_metadata_cluster %>%
      filter(cluster == 1) %>%
      mutate(cluster2 = DRAM_products_metadata_cluster_clust) %>%
      group_by(cluster2) %>%
      slice_sample(n = 1) %>%
      select(genome,cluster2)

#Add whether selected to metadata
DRAM_products_metadata_cluster2 <- DRAM_products_metadata_cluster %>%
        left_join(DRAM_products_metadata_cluster_selected,by="genome") %>%
        relocate(cluster,cluster2) %>%
        mutate(cluster2 = ifelse(is.na(cluster2), 0, 1))

write.table(DRAM_products_metadata_cluster_selected,"data/selected_genomes.tsv",row.names=F)

#Color PCA
DRAM_products_pc %>%
  augment(DRAM_products_metadata_cluster2) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color=as.factor(cluster2))) +
    geom_point(size = 0.5, alpha=0.5) +
    scale_color_manual(values=c("#cccccc", "#cc14c6"))+
    theme_minimal() +
    theme(legend.position = "none")

#Color tsne
DRAM_products_tsne$Y %>%
  as_tibble() %>%
  mutate(genome=DRAM_products_metadata_cluster2$genome) %>%
  left_join(DRAM_products_metadata_cluster2, by="genome") %>%
  rename(tSNE1="V1", tSNE2="V2") %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color=as.factor(cluster2)))+
    geom_point(size=0.5, shape=16, alpha=0.8) +
    scale_color_manual(values=c("#cccccc", "#cc14c6"))+
    theme_minimal() +
    theme(legend.position = "none")

#####
# Variable distribution
#####

#All genomes
DRAM_products_metadata_cluster2 %>%
      select(coding_density,gc_percentage,genome_size,protein_count) %>%
      gather() %>%
      ggplot(., aes(value)) +
        geom_histogram(bins = 30) +
        facet_wrap(~key, scales = 'free_x')

#Selected genomes
DRAM_products_metadata_cluster2 %>%
      filter(cluster2 == 1) %>%
      select(coding_density,gc_percentage,genome_size,protein_count) %>%
      gather() %>%
      ggplot(., aes(value)) +
        geom_histogram(bins = 30) +
        facet_wrap(~key, scales = 'free_x')
