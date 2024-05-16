---
title: "SA Phage Rotation Project"
author: "Sophia Aredas"
date: "15 May, 2024"
output:
  html_document:
    code_folding: show
    highlight: default
    keep_md: yes
    theme: journal
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: console
---

```{=html}
<style>
pre code, pre, code {
  white-space: pre !important;
  overflow-x: scroll !important;
  word-break: keep-all !important;
  word-wrap: initial !important;
}
</style>
```


```r
# For width of code chunks and scroll bar
knitr::opts_chunk$set(eval = TRUE,
                      echo = TRUE,
                      include = TRUE,
                      warning = FALSE,
                      collapse = FALSE,
                      message = FALSE,
                      dpi=200, dev = "png",
                      error = TRUE,
                      fig.path="figures/",
                      fig.align = "center")
```


```r
# Efficiently load packages
pacman::p_load(readr, stringr, tidyr, tidyverse, dplyr, ggplot2, tidytext, scales, wesanderson, ggpubr, wacolors, rcartocolor, treeio, ggtreeExtra, phyloseq, ggtree, ggstar, phytools, install = FALSE)
```

### Load summary files produced by geNomad


```r
#set working directory
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")

#1. looking at plasmid genes
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
plasmid_genes <- read.table(file = "plasmid_genes.tsv", sep ="\t", header = TRUE)

#2. looking at plasmid summary 
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
plasmid_summary <- read.table(file = "plasmid_summary.tsv", sep ="\t", header = TRUE)

#3 looking at virus genes
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
virus_genes <- read.table(file = "virus_genes.tsv", sep ="\t", header = TRUE)

#4 looking at virus summary
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
virus_summary <- read.table(file = "virus_summary.tsv", sep ="\t", header = TRUE)
```

## Creating data frames so that we can analyze the geNomad outputs

### Creating dataframe of plasmid genes


```r
#1. plasmid genes
#need to create dataframe
df_plasmidgenes <- data.frame(plasmid_genes)


#split the first column called gene into multiple columns
#because there is so much going on with the column name we are going to separate each word into its own column and then mutate it together (dont be alarmed)
df_plasmidgenes <- df_plasmidgenes %>% separate(gene, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'plasmid', 'genes', 'tsv', 'Contig', 'bp_length', 'Contig_Number'))

#now we are going to mutate the combine the columns together
df_plasmidgenes$Seq_Name <- paste(df_plasmidgenes$Seq_Name, df_plasmidgenes$DASTool, df_plasmidgenes$Bin, sep = "_")
df_plasmidgenes$fa_info <- paste(df_plasmidgenes$fa_info, df_plasmidgenes$plasmid, df_plasmidgenes$genes, df_plasmidgenes$tsv, sep = "_")
df_plasmidgenes$Contig <- paste(df_plasmidgenes$Contig, df_plasmidgenes$bp_length, sep = "_")


#then we need to get rid of the extra columns that we do not need
#df_plasmidgenes <- subset(df_plasmidgenes, select = -c(DASTool, plasmid, genes, tsv, bp_length)) #invalid argument need to fix

#view to make sure it looks good
```

### Creating dataframe of plasmid summary


```r
#2. plasmid summary
#need to create dataframe
df_plasmidsum <- data.frame(plasmid_summary)

#split the first column called gene into multiple columns
df_plasmidsum <- df_plasmidsum %>% separate(seq_name, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'plasmid', 'genes', 'tsv', 'Contig', 'bp_length'))

#now we are going to mutate the combine the columns together
df_plasmidsum$Seq_Name <- paste(df_plasmidsum$Seq_Name, df_plasmidsum$DASTool, df_plasmidsum$Bin, sep = "_")
df_plasmidsum$fa_info <- paste(df_plasmidsum$fa_info, df_plasmidsum$plasmid, df_plasmidsum$genes, df_plasmidsum$tsv, sep = "_")
df_plasmidsum$Contig <- paste(df_plasmidsum$Contig, df_plasmidsum$bp_length, sep = "_")

#then we need to get rid of the extra columns that we do not need
df_plasmidsum <- subset(df_plasmidsum, select = -c(DASTool, plasmid, genes, tsv, bp_length)) 

#view to make sure it looks good
```

### Creating virus genes dataframe


```r
#3. virus genes
#need to create dataframe
df_virusgenes <- data.frame(virus_genes)

#split the first column called gene into multiple columns
df_virusgenes <- df_virusgenes %>% separate(gene, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'virus', 'genes', 'tsv', 'Contig', 'bp_length', 'Contig_Number'))

#now we are going to mutate the combine the columns together
df_virusgenes$Seq_Name <- paste(df_virusgenes$Seq_Name, df_virusgenes$DASTool, df_virusgenes$Bin, sep = "_")
df_virusgenes$fa_info <- paste(df_virusgenes$fa_info, df_virusgenes$virus, df_virusgenes$genes, df_virusgenes$tsv, sep = "_")
df_virusgenes$Contig <- paste(df_virusgenes$Contig, df_virusgenes$bp_length, sep = "_")

#then we need to get rid of the extra columns that we do not need
df_virusgenes <- subset(df_virusgenes, select = -c(DASTool, virus, genes, tsv, bp_length)) #invalid argument need to fix

#view to make sure it looks good
```

### Creating virus summary dataframe


```r
#4. virus summary
#need to create dataframe
df_virussum <- data.frame(virus_summary)

#split the first column called gene into multiple columns
df_virussum <- df_virussum %>% separate(seq_name, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'virus', 'genes', 'tsv', 'Contig', 'bp_length'))

#now we are going to mutate the combine the columns together
df_virussum$Seq_Name <- paste(df_virussum$Seq_Name, df_virussum$DASTool, df_virussum$Bin, sep = "_")
df_virussum$fa_info <- paste(df_virussum$fa_info, df_virussum$virus, df_virussum$genes, df_virussum$tsv, sep = "_")
df_virussum$Contig <- paste(df_virussum$Contig, df_virussum$bp_length, sep = "_")

#then we need to get rid of the extra columns that we do not need
df_virussum <- subset(df_virussum, select = -c(DASTool, virus, genes, tsv, bp_length)) 

#view to make sure it looks good
```

## Combining Gus's data frame and mine together to map the bin to the fraction class

### The names on the 346 MAG files are from the filtration steps. The free coassembly was done from the free filtration, particle coassembly was done from the particle, and DASTools was presumably done from the whole filtration


```r
#but the names on each of the files are meaningless. The file could be like free-coassembly bin 90 but we do not know if they are actually FL, PA, or generalists.
#Thanks to Gus he provided MAG_Growth_Fraction.RData which there were calculations to see if bins prefer to be FL, PA, or Generalists. I will be adding his column into my dataset so that way I am able to compare between bacterial lifestyles

#load in the MAG_Growth_Fraction.RData to your environment
load("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/MAG_Growth_and_Fraction.RData")

# Now lets read in our checkm file which has bin lengths. We need this order to correctly calculate genes per bin
checkm <- read.csv("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/summarized_coverage_perbin.csv")

checkm_binlength <- checkm %>% 
  group_by(bin, bin_length) %>% 
  summarise(.groups = "drop")

#left join to add to mag_full_df 
mag_full_df <- left_join(mag_full_df, checkm_binlength, by = c("bin" = "bin"))

#gsub so that free and particle coassemblies match
mag_full_df$bin <- gsub("-", "_", mag_full_df$bin)

#now we also need to read in the checkM_bins file so that we can get our completeness score
checkm_completeness <- read.table(file = "/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/checkM_bins826_50completeness.tsv", sep ="\t", header = TRUE)

#gsub so that free and particle coassemblies match
checkm_completeness$bin_ID <- gsub("-", "_", checkm_completeness$bin_ID)
#now delimit to get rid of ".contigs"
checkm_completeness <- data.frame(checkm_completeness)

#split the first column called gene into multiple columns
checkm_completeness <- checkm_completeness %>% separate(bin_ID, c('Seq_Name', 'DASTool', 'Bin', 'Contigs'))
#then we will take the multiple columns and paste them togther to get rid of the .contigs
checkm_completeness$Seq_Name <- paste(checkm_completeness$Seq_Name, checkm_completeness$DASTool, checkm_completeness$Bin, sep = "_")
#then remove the unnecessary columns
checkm_completeness <- subset(checkm_completeness, select = -c(DASTool, Bin, Contigs)) 
#now we are just selecting the columns that we want
checkm_completeness <- checkm_completeness %>% 
  select(c(Seq_Name, Completeness))

#now lets join the checkm_binglength with the checkm_completeness to mag_full_df that already has the bin_length
mag_full_df <- inner_join(checkm_completeness, mag_full_df, by = c("Seq_Name" = "bin"))

#the bin_length column was created by summing all the lengths of the contigs within the bin together 
#now lets go from bp to Mbp
mag_full_df <- mag_full_df %>% 
  mutate(length_mb = bin_length / 1000000)

#now to get the "true" length of the genome we will multiple the length_mb and multiply it by the completeness score
#the estimated length is now in Mb
mag_full_df <- mag_full_df %>% 
  mutate(estimated_length = length_mb * Completeness)

#join columns together with merge left_join
merge_virussum_magfull <- left_join(df_virussum, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
merge_plasmidgenes_magfull <- left_join(df_plasmidgenes, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
merge_virusgenes_magfull <- left_join(df_virusgenes, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
merge_plasmidsum_magfull <- left_join(df_plasmidsum, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
```

## Time to visualize

#### Relative Abundance of Fraction Class

##### Bar Graph with no percentages just the relative abundance count

```r
#i am making a summary df to calculate the sum of fraction_class
#looking at Gus's MAG_full_df
#this is the bargraph version with no percentages
count_magfull <- mag_full_df %>% 
  group_by(fraction_class) %>% 
  summarise(count = n(), .groups = "drop")

count_fraction_class <- count_magfull %>% 
  ggplot(aes(x = fraction_class, y = count, fill = fraction_class))+
  geom_bar(stat = "identity", position = "dodge")
  
count_fraction_class
```

<img src="figures/relative abundance of fraction_class-1.png" style="display: block; margin: auto;" />

##### Bar Graph with percentages


```r
#this is the percentage bar graph version
percent_mag_full <- mag_full_df %>% 
  group_by(fraction_class) %>% 
  summarise(count = n()) %>% 
  mutate(percent = count/sum(count))


# ggplot(percent_mag_full, aes(x = fraction_class, y = percent, fill = fraction_class, label = scales::percent(percent(stat(count))))) +
#   geom_bar(stat = "identity") +
#   geom_text(stat = "identity",
#             position = position_dodge(.9),
#             vjust = -0.5,
#             size = 3) +
#   labs(x = "Lifestyle", y = "Percent", fill = "Fraction Class") +
#   theme_minimal(base_size = 14)
# 
# percent <- ggplot(percent_mag_full, aes(x = fraction_class, y = percent, fill = fraction_class))+
#   geom_bar(stat = "identity", position = "dodge")
# 
# percent +
#   scale_y_continuous(labels = scales::percent,
#                      breaks= scales::pretty_breaks(n = 7))

#then all together
percent <- ggplot(percent_mag_full, aes(x = fraction_class, y = percent, fill = fraction_class, label = scales::percent(percent)))+
  geom_col(width = 0.5)+
  labs(x = "Lifestyle",
       y = "Percentage of MAGs")+
  scale_y_continuous(labels = scales::percent,
                     breaks= scales::pretty_breaks(n = 7))+
  geom_text(nudge_y = -.03,
            color = "white",
            size = 5, 
            fontface = "bold")+
  scale_fill_wa_d("rainier", reverse=TRUE)+
  theme_classic()
percent
```

<img src="figures/percentage bar graph version-1.png" style="display: block; margin: auto;" />

```r
#percent + labs(fill = "Lifestyle")

#009392,#39b185,#9ccb86,#e9e29c,#eeb479,#e88471,#cf597e# A tibble: 3 × 2
 # fraction_class count
 # <fct>          <int>
#1 Free              63
#2 Particle          52
#3 Generalist       231
```
### Plasmid Analysis

#### Lets look at the breakdown of plasmids in our data set

```r
plasmid_presence <- merge_plasmidsum_magfull

#now lets normalize
plasmid_presence <- plasmid_presence %>%
  group_by(Seq_Name, fraction_class, estimated_length) %>%
  summarise(
    count_plasmid_hits = n()
  )

#now we can divide number of hits by estimated length
plasmid_presence$abundance <- plasmid_presence$count_plasmid_hits / plasmid_presence$estimated_length


#now lets visualize the relative abundance I have not normalized it to MAG
plas_presence <- ggplot(plasmid_presence, aes(x = fraction_class, y= abundance, fill = fraction_class))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point()+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE)+
    stat_compare_means(comparisons = list(
      c("Free", "Particle"),
      c("Free", "Generalist"),
      c("Particle", "Generalist")),
                       method = "wilcox.test", paired = FALSE, label = "p.adj.format",hide.ns = TRUE,
                       p.adjust.method = "fdr") +
  xlab("Lifestyle")+
  ylab("# Plasmid Hits / MAG")+
  theme_classic()

plas_presence
```

<img src="figures/plasmid breakdown-1.png" style="display: block; margin: auto;" />

##### Visualize the number of plasmid encoded genes and fraction class


```r
#looking at the number of peg divided by mag but represented in a box plot
#so we will go about this by dividing the number of genes column by the length of bin column
plasmid_sum_magfull_gene_count <- merge_plasmidsum_magfull%>% 
  group_by(Seq_Name, Bin, fraction_class, Phylum, estimated_length) %>%
  summarise(count_n_genes = sum(n_genes))


#now we will divide the number of genes per length to get the abundance per assembly
plasmid_sum_magfull_gene_count$genes_per_bin <- plasmid_sum_magfull_gene_count$count_n_genes / plasmid_sum_magfull_gene_count$estimated_length

#now lets visualize our data
peg_phy <- ggplot(plasmid_sum_magfull_gene_count, aes(x = fraction_class, y = genes_per_bin, fill = Phylum))+
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_carto_d(palette = "Geyser")+
  xlab("Lifestyle")+
  ylab("Number of Plasmid Genes Encoded per Mbp per MAG ")+
  theme_classic()

peg_phy
```

<img src="figures/visualize plasmid n_genes and fraction_class-1.png" style="display: block; margin: auto;" />

```r
#now lets look at plasmid encoded genes per bin 
peg_per_bin <- ggplot(plasmid_sum_magfull_gene_count, aes(x = fraction_class, y = genes_per_bin, fill = fraction_class))+
  geom_boxplot(outlier.shape = NA, alpha = .8)+
  geom_point()+
  xlab("Lifestyle")+
  ylab("Number of Plasmid Genes Encoded per Mbp per MAG")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE)+
  stat_compare_means(comparisons = list(
      c("Free", "Particle"),
      c("Free", "Generalist"),
      c("Particle", "Generalist")),
                       method = "wilcox.test", paired = FALSE, label = "p.adj.label",hide.ns = FALSE,
                       p.adjust.method = "fdr") +
  theme_classic()

peg_per_bin
```

<img src="figures/visualize plasmid n_genes and fraction_class-2.png" style="display: block; margin: auto;" />

#### Density plot looking at number of plasmid encoded genes/MAG (x axis) and the density (y axis)


```r
peg_density <- ggplot(plasmid_sum_magfull_gene_count, aes(x = genes_per_bin, color = fraction_class, fill = fraction_class))+
  geom_density(alpha = 0.3)+
  xlab("Number of Plasmid Genes Encoded per Mbp per MAG")+
  ylab("Density")+
  scale_fill_manual(
    values = c("#DF3383","#759C44","#9FB6DA")) +
  theme_classic()

peg_density
```

<img src="figures/density plot of plasmid encoded genes (number of peg per mag)-1.png" style="display: block; margin: auto;" />

#### Visualizing conjugation genes divided by MAG


```r
#now we will visualize the amount of conjugation genes within the merge_plasmidsum dataframe
congenee <- merge_plasmidgenes_magfull

#now we will make all conjugation type genes as an umbrella
congenee <- congenee %>% 
  filter(!is.na(annotation_conjscan)) %>%  
  mutate(conj_gene_type = case_when(grepl("trb", annotation_conjscan) ~ "trb genes",
                                    grepl("vir", annotation_conjscan) ~ "vir genes",
                                    grepl("MOB", annotation_conjscan) ~ "MOB genes",
                                    grepl("tra", annotation_conjscan) ~ "tra genes",
                                    grepl("t4cp1", annotation_conjscan) ~ "tra genes")) %>% 
  select(Seq_Name, Bin, annotation_conjscan, conj_gene_type, fraction_class, estimated_length, Phylum, log2FoldChange)

# now we can visualize the number of conjugation genes per bin
# count all the genes in the bin regardless of conjugation gene type
congenee_all <- congenee %>%
  group_by(Seq_Name, fraction_class, estimated_length, log2FoldChange) %>%
  summarise(count_conjugation = n())
#this is counting bin but taking into account conjugation gene type
congenee_type <- congenee %>%
  group_by(Seq_Name, conj_gene_type, fraction_class, Phylum, estimated_length, log2FoldChange) %>%
  summarise(
    Count_all_conj_genes_bin = n()
  )

#now lets get the relative abundance of all the conjugation genes that are not specified by conjugation gene type
congenee_all$conj_genes_per_bin <- congenee_all$count_conjugation / congenee_all$estimated_length

#but if we want to see the amounts of conjugation gene type per bin and unique type of conjugation gene 
congenee_type$conj_genes_per_bin <- congenee_type$Count_all_conj_genes_bin/ congenee_type$estimated_length

#now lets try to merge the conj_genes_per_bin to get the conugation gene type
congenee_tgt <- left_join(congenee_all, congenee_type, by = c("Seq_Name" = "Seq_Name"))
#now that we have ourgenes/sequence name we can plot our data to see what is happening
conj_genes_per_bin <- ggplot(congenee_all, aes(x = fraction_class, y = conj_genes_per_bin, fill = fraction_class))+
  geom_boxplot() +
  geom_point()+
  xlab("Lifestyle")+
  ylab("Number of Conjugation Genes per Mbp per MAG")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE)+
  # stat_compare_means(comparisons = list(
  #     c("Free", "Particle"),
  #     c("Free", "Generalist"),
  #     c("Particle", "Generalist")),
  #                      method = "wilcox.test", paired = FALSE, label = "p.signif",hide.ns = FALSE,
  #                      p.adjust.method = "fdr") +
  theme_classic()
conj_genes_per_bin
```

<img src="figures/visualize conjugation genes divided by mag-1.png" style="display: block; margin: auto;" />

```r
#now lets visualize the type of conjugation gene in each of the bins
#warning some bins have multiple kinds of conjugation genes so it may look weird but thats ok hopefully 
conj_genes_per_bin_type <- ggplot(congenee_type, aes(x = fraction_class, y = conj_genes_per_bin, fill = conj_gene_type))+
  geom_bar(position = "dodge", stat = "identity") +
  xlab(" Bacterial Lifestyle")+
  ylab("Number of Conjugation Genes per Mbp per MAG")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE) +
  labs(fill = "Conjugation Genes") +
  theme_classic()
conj_genes_per_bin_type
```

<img src="figures/visualize conjugation genes divided by mag-2.png" style="display: block; margin: auto;" />

```r
# plasmid_sum_magfull_gene_count %>%
#   group_by(fraction_class, genes_per_bin) %>%
#   get_summary_stats(genes_per_bin)

#now lets use ggsave to save this plot at a higher dpi and size
ggsave("conj_genes_per_bin_type_plot.png", plot = conj_genes_per_bin_type, width = 5, height = 4, dpi = 1000)


#I am going out of order but here I am making a plot to look at the differential abundance of MAGs, phyla, and lifestyle to see where conjugations genes are distributed
log2fold_mags <- ggplot(mag_full_df, aes(x = log2FoldChange, y = Phylum, fill = fraction_class, shape = fraction_class)) +
  geom_col(position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_point(alpha = 0.5, position = position(width = 0.4), show.legend = TRUE) +
  scale_fill_manual(values = fraction_class_colors) +
  scale_shape_manual(values = c("Free" = 24, "Generalist" = 21, "Particle" = 25)) +
  theme_classic() +
  theme(legend.position = "bottom") 
```

```
## Error in position(width = 0.4): could not find function "position"
```

```r
log2fold_mag
```

```
## Error in eval(expr, envir, enclos): object 'log2fold_mag' not found
```

```r
log2fold_conjugation <- ggplot(congenee_type, aes(x = log2FoldChange, y = Phylum, fill = conj_gene_type, shape = fraction_class)) +
  geom_col(position = position_dodge(width = .95), show.legend = FALSE) +
  geom_point(alpha = 0.5, position = position_jitterdodge(), show.legend = TRUE) +
  xlab("Log2FoldChange")+
  ylab("Phylum")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE) +
  scale_shape_manual(values = c("Free" = 24, "Generalist" = 21, "Particle" = 25)) +
  theme_classic() 
log2fold_conjugation
```

<img src="figures/visualize conjugation genes divided by mag-3.png" style="display: block; margin: auto;" />

```r
#this looks insane I need to talk to marian tmr

#ggsave("log2foldmag.png", plot = log2fold_mags, width = 5, height = 4, dpi = 1000)
```

### Phage Analysis

#### Lets first do a bar graph visualizing of viruses + unclassified  + NA in sample

```r
tax_phass <- merge_virussum_magfull
#split the column called taxonomy.filename into multiple columns by the delimiter ";"
tax_phass <- tax_phass %>% separate(taxonomy.filename, c('Vir_Domain', "Realm", "Kingdom", "Phylum_vir", "Class_cvir", "Order_vir", "Species_vir"))

#now lets normalize
phage_presence <- tax_phass %>%
  group_by(Seq_Name, fraction_class, estimated_length, Vir_Domain, Phylum_vir, Phylum) %>%
  summarise(
    count_virus_hits = sum(n_genes)
  )

#but lets make the NA unclassified ones into unclassified
phage_presence$Phylum_vir[is.na(phage_presence$Phylum_vir)] <- "Unclassified" 

#now we can divide number of hits by estimated length
phage_presence$abundance <- phage_presence$count_virus_hits / phage_presence$estimated_length


#now lets visualize the relative abundance thta is normalized for number of phage hits divided by mag
virus_presence <- ggplot(phage_presence, aes(x = fraction_class, y= abundance, fill = fraction_class))+
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_point()+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE)+
    stat_compare_means(comparisons = list(
      c("Free", "Particle"),
      c("Free", "Generalist"),
      c("Particle", "Generalist")),
                       method = "wilcox.test", paired = FALSE, label = "p.adj.format",hide.ns = TRUE,
                       p.adjust.method = "fdr") +
  xlab("Lifestyle")+
  ylab("Number of Phage-Encoded Genes per Mbp per MAG")+
  theme_classic()

virus_presence
```

<img src="figures/virus presence-1.png" style="display: block; margin: auto;" />


#### Bar graph visualizing the number of phage genes and phage taxonomy


```r
#we will delmit so that we can see what phages are present in our dataset 
#create different variable of the df so that in case we mess up its gonna be ok
tax_phass <- merge_virussum_magfull

#split the column called taxonomy.filename into multiple columns by the delimiter ";"
tax_phass <- tax_phass %>% separate(taxonomy.filename, c('Vir_Domain', "Realm", "Kingdom", "Phylum_vir", "Class_cvir", "Order_vir", "Species_vir"))

#changing in the phylum_vir column the names "NA" to unclassified
tax_phass$Phylum_vir[is.na(tax_phass$Phylum_vir)] <- "Unclassified"

#lets combine with magfull_df
tax_phass <- left_join(tax_phass, mag_full_df, by=c("Seq_Name"="Seq_Name"))
#lets use merge to combine tax_phass and count_phab_seq
phab_bac_tax <- full_join(phage_presence, mag_full_df, by = c("Seq_Name"= "Seq_Name"))

#now we can visualize the number of phage genes and look at the phage taxonomy
number_of_phage_genes <- ggplot(phage_presence, aes(x = fraction_class, y = abundance, fill = Phylum_vir))+
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_carto_d(palette = "Geyser")+
  #scale_fill_wa_d(palette = "rainier", reverse = FALSE)+
  xlab("Lifestyle")+
  ylab("Number of Phage Genes Encoded / MAG")+
  theme_classic()

number_of_phage_genes 
```

<img src="figures/looking at the phage taxonomy that is on row and it is separated by ";"-1.png" style="display: block; margin: auto;" />


#### Number of phage encoded genes per MAG


```r
#now lets divide the number of genes by MAG (AKA lest visualize the number of phage
#count all the genes in the bin 
phab_genes <- merge_virussum_magfull
phab_genes <- phab_genes %>%
  group_by(Seq_Name, Bin, fraction_class, estimated_length, Phylum) %>% 
  summarise(
    .groups = "drop",
    genes_per_bin_phage = sum(n_genes))

#now we will divide the number of genes per length to get the abundance per assembly
phab_genes$genes_per_bin_phage <- phab_genes$genes_per_bin_phage / phab_genes$estimated_length


#now lets plot to see amount of phage encoded genes per bin by fraction class
phab_gene_per_bin <- ggplot(phab_genes, aes(x = fraction_class, y = genes_per_bin_phage, fill = fraction_class))+
  geom_boxplot(outlier.shape = NA, alpha = .8)+
  geom_point()+
  xlab("Lifestyle")+
  ylab("Number of Phage-Encoded Genes per Mbp per MAG")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE)+
  stat_compare_means(comparisons = list(
      c("Free", "Particle"),
      c("Free", "Generalist"),
      c("Particle", "Generalist")),
                       method = "wilcox.test", paired = FALSE, label = "p.adj.format", hide.ns = TRUE,
                       p.adjust.method = "fdr") +
  theme_classic()
phab_gene_per_bin
```

<img src="figures/phage genes divided per MAG-1.png" style="display: block; margin: auto;" />

#### Density plot of phage encoded genes per MAG


```r
phage_encoded_density <- ggplot(phab_genes, aes(x = genes_per_bin_phage, color = fraction_class, fill = fraction_class))+
  geom_density(alpha = 0.3)+
  xlab("Number of Phage Encoded Genes per Mbp per MAG")+
  ylab("Density")+
  scale_fill_manual(
    values = c("#DF3383","#759C44","#9FB6DA")) +
  theme_classic()

phage_encoded_density
```

<img src="figures/Density plot of phage encoded genes per MAG-1.png" style="display: block; margin: auto;" />

#### Looking at bacterial lifestyle (x axis) and the number of phage encoded genes / MAG faceted by phage phylum/ taxonomy


```r
#changing in the phylum_vir column the names "NA" to unclassified
phab_bac_tax$Phylum_vir[is.na(phab_bac_tax$Phylum_vir)] <- "Unclassified"

#now plotting using a facet to see from the most common phages and where they are present in bacterial phyla 
#ommitted phages that were present in only one bacterial phylum
number_of_phage_genes1 <- 
  phage_presence %>%
  dplyr::filter(Phylum_vir %in% c("Uroviricota", "Nucleocytoviricota", "Unclassified")) %>%  
  ggplot(aes(x = Phylum, y = abundance, fill = Phylum))+
  facet_wrap(~ Phylum_vir) + 
  geom_bar(position = "dodge", stat = "identity") +
 theme_classic() + 
   theme(
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(2, "lines")
    )+
  scale_fill_carto_d(palette = "Geyser")+
  labs(x = "Lifestyle", 
       y = "Number of Phage Genes Encoded per Mbp per MAG");number_of_phage_genes1
```

<img src="figures/faceting virus phylum and bacteria phylum to see which phages are found in bacteria-1.png" style="display: block; margin: auto;" />

```r
 #now looking at lifestyle as opposed to bacterial phylum
number_of_phage_genes2 <- 
  phage_presence %>%
  #dplyr::filter(Phylum_vir %in% c("Uroviricota", "Nucleocytoviricota", "Unclassified")) %>%  
  ggplot(aes(x = Phylum_vir, y = abundance, fill = Phylum_vir))+
  facet_wrap(~ fraction_class) + 
  geom_bar(position = "dodge", stat = "identity") +
 theme_classic() + 
   theme(
    #axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(2, "lines")
    )+
  #scale_fill_wa_d(palette = "rainier", reverse = FALSE)+
  scale_fill_carto_d(palette = "Geyser")+
  labs(x = "Lifestyle", 
       y = "Number of Phage Genes Encoded per Mbp per MAG");number_of_phage_genes2
```

<img src="figures/faceting virus phylum and bacteria phylum to see which phages are found in bacteria-2.png" style="display: block; margin: auto;" />

### Using GTDBTK to do phylogenetic analysis

#### The code that I used is located in my code directory called gtdbtk.sh


```gtdbtk
#To use gtdbtk there is no need to install it as it is a huge download and we already have it on our server located in databases. Please go to /local/workdir/databases/gtdbtk_release214/release214/release214/taxonomy to use this
source $HOME/miniconda3/bin/activate

#Set GTDBTK_DATA_PATH to new reference db directory, for example
#if you do gtdbtk-2.1.1 that is the release version for the release207_v2
#if you want to do the most updated release (to me that seems more reasonable) then you will want to use the gtdbtk-2.3.

conda create -n gtdbtk-2.3.2 -c conda-forge -c bioconda gtdbtk=2.3.2

# To activate this environment, use
#
#     $ conda activate gtdbtk-2.3.2
#
# To deactivate an active environment, use
#
#     $ conda deactivate

#now activate the conda environment so that we can add the stuff to the environment 
conda activate gtdbtk-2.3.2
#The conda package is bundled with a script download-db.sh (source) that will automatically download, and extract the GTDB-Tk reference data. The script will be on the system path so simply run:
#this part takes FOREVER just fyi use screen or something
download-db.sh

#set path 
conda env config vars set GTDBTK_DATA_PATH="/home/sna49/miniconda3/envs/gtdbtk-2.3.2/share/gtdbtk-2.3.2/db"

#checking which version just in case (keep in mind this part takes a minute)
gtdbtk check_install

#created new directory called gtdbtk 
/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/346_MAGs
  #I recopied the original 346 MAGs but jsut the .fa version
  
#we will now begin to do our analysis with the command called classify_wf
  #this workflow consists fo 4 steps:
    #ani_screen: compares our genomes against a MASH database composed of all GTDB representative genoems then verifies the best mash hits using FASTANI which is a method used to estimate average nucleotide identity using alignment free approximate sequence mapping. it is good for finished and draft genomes
    #identify: uses prodigal and uses HMM mdoels and HMMER package to identify the 120 bacterial and 53 archaeal marker genes used for phylogenetic inference. MSA obtained by aligning amrker gene to respective HMM model
    #align: concatenates the aligned marker genes and filters the concatenated MSP to approx 5k amino acids
    #classify: uses pplacer to find the ML placement of each genome in teh gtdbtk reference tree. GTDBTK classifies each genome based on its placment in refernece tree, relative evolutionary divergence, and/or ANI to reference geomes
    
    #workflow command: gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
      #genomes must be in FASTA format but .gz is acceptable which is awesome bc thats what we got
      
gtdbtk classify_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/dereplicated_bins/346_MAGs --out_dir classify_gtdbtk_output_skip --skip_ani_screen -x .fa --cpus 35
#gtdbtk = command
#classify_wf = pipeline used
##--genome_dir = pointing to the genome directory or wherever our files are at
#out_dir = pointing to the output directory
#classify_gtdbtk_output = this is the output file that is being created
#--mash_db 2.3.2 = this version of GTDBTK requires that you input the mash_db version (2.3.2 is the one we downloaded). This Using the --mash_db option will indicate to GTDB-Tk the path of the sketched Mash database require for ANI screening.
  #If no database are available ( i.e. this is the first time running classify ), the --mash_db option will sketch a new Mash database that can be used for subsequent calls.
#-x .fa = -x means extension files. the files we have are .fa
#--cpus 100 = these are the cpus that are being used by the server to calculate

#general notes about downloading
  #does not seem to like the .gz files 
  #instead I just did the original .fa files 

  
#so far GTDBTK gives us the msa alignment when we use the code above, this can be used for an unrooted tree which is totally fine

#but we can also try the de novo wf to give us a rooted tree and we would use planctomycetes as our outgroup
#i didnt end up using the rooted tree i did not have the capacity to think about this lol
gtdbtk de_novo_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/346_MAGs --outgroup_taxon p__Chloroflexota --bacteria  --out_dir de_novo_wf_test1 -x .fa --cpus 100
```

#### Installing fastree


```bash_fasttree
#There are 2 previous version(s) available. Any version can be accessed either by typing full path, or by adding it to the PATH and then typing its name at the prompt. NOTE: you need to set PATH only once per login/session.

#add to path
export PATH=/programs/FastTree-2.1.11:$PATH

#then type name to get the program activated
FastTree

# http://www.microbesonline.org/fasttree/#Usage link for running the code 

#when looking at our results we see that the files are in AA format
chmod +rwx Data/gtdbtk/classify_gtdbtk_output/align
chmod +rwx  data/MAG_gtdb_magonly_tree/align/gtdbtk.bac120.user_msa.fasta.gz
gunzip data/MAG_gtdb_magonly_tree/align/gtdbtk.bac120.user_msa.fasta.gz
FastTree data/MAG_gtdb_magonly_tree/align/gtdbtk.bac120.user_msa.fasta > data/MAG_gtdb_magonly_tree/gtdb_magonly.tree

#for a protein alignment #FastTree input_aln.fasta > input.tre FOR UNROOTED TREE
chmod +rwx /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output/align/gtdbtk.bac120.user_msa.fasta.gz
gunzip /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output/align/gtdbtk.bac120.user_msa.fasta.gz
FastTree /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output/align/gtdbtk.bac120.user_msa.fasta > /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output/tree.tree
conda deactivate

#this code is for our rooted tree created by de_novo_wf
chmod +rwx /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output_skip/align/gtdbtk.bac120.user_msa.fasta.gz
gunzip /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output_skip/align/gtdbtk.bac120.user_msa.fasta.gz
FastTree /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output_skip/align/gtdbtk.bac120.user_msa.fasta > /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output_skip/align/unrooted.tree
```

#### Creating phyloseq object


```r
#to create phyloseq object we need 3 things: OTU, Taxonomy, Samples
#I am using Gus's phyloseq object and adding my own taxonomy to it as GTDBTK was able to update the taxonomy (Thank you Gus!!)

#brining in the relative abundance table. Gus calculated with bwa, bbmap:pileup, and "summarizing_percontig_coverage.R"
cov_table_sum<-read.csv ("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/phyloseq/transposed_norm_avg_coverage_perbin.csv")
#use gsub clean up data NEED TO DO THIS PART I LEFT OFF HERE
cov_table_sum$bin <- gsub("-", "_", cov_table_sum$bin)

#define the row names from OTU column 
matrix_cov_table_sum<-as.matrix.data.frame(cov_table_sum[3:18])
rownames(matrix_cov_table_sum)<-cov_table_sum$bin



#then bring in the sample matrix
metadata<-read.csv("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/phyloseq/sample_metadata.csv")%>%
  mutate(sample = paste0("X", sample_ID))%>%
  select(-sample_ID)

rownames(metadata)<-metadata$sample

#bring in our taxonomy table 
# need to clean up the data as well
taxadata <- read.delim("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/classify_gtdbtk_output_skip/gtdbtk.bac120.summary.tsv") %>% 
  select(user_genome, classification) %>% 
  separate(classification, sep = ";", into = c("Domain","Phylum","Class","Order","Family","Genus","Species"))%>%
  mutate(across(!user_genome, function(x)sub(x, pattern = ".__", replacement = "")),
         user_genome = stringr::str_remove(user_genome, ".contigs"))

#use gsub clean up data 
taxadata$user_genome <- gsub("-", "_", taxadata$user_genome)

#adding in our mag_full_df file
mag_full <- mag_full_df 

mag_full <- left_join(mag_full, taxadata, by = c("Seq_Name" = "user_genome"))
#cleaning up our merged file
mag_full <- mag_full %>% 
  select(-c("Domain.x", "Phylum.x", "Class.x", "Order.x", "Family.x", "Genus.x", "Species.x", "CUBHE", "GC", "GCdiv", "ConsistencyHE", "CUB", "CPB", "nHE", "dCUB", "d", "LowerCI", "UpperCI", "FilteredSequences", "log2FoldChange", "padj", "growth_class",))

#Bring in gRodon data
growthrate<-read.csv("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/phyloseq/bakta_mags_growthrates.csv")%>%
  select(2:Sample)%>%
  mutate(Sample = str_remove(Sample, pattern = ".contigs"))

growthrate$Sample <- gsub("-", "_", growthrate$Sample)

taxadata<-taxadata%>%
  inner_join(growthrate, by = c("user_genome" = "Sample"))


taxa_matrix<-as.matrix(taxadata[2:ncol(taxadata)])
rownames(taxa_matrix)<-taxadata$user_genome

#Bringing in our phylogenetic tree. this is the one that is created from our classify_wf with the --skip_ani_screen flag
tree<-read.newick("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/phyloseq/unrooted.tree")%>%
  as_tibble()%>%
  mutate(label = stringr::str_remove(label, ".contigs"))%>%
  as.phylo()

tree$tip.label <- gsub("-", "_", tree$tip.label)

#Construct our phyloseq sample, otu, and taxa tables
SAMP<-sample_data(metadata)
OTU_sum<-otu_table(matrix_cov_table_sum, taxa_are_rows = T)
TAX<-tax_table(taxa_matrix)

#Put them together into a phyloseq object
mag_phylo_sa <- phyloseq(OTU_sum, SAMP, TAX, tree)

#Saving to RData file for analysis
save(mag_phylo_sa, file = ("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/MAG_phyloseq_SA.RData"))
```

#### Using ggtree and other packages to view unrooted phylogenetic tree

```r
#read in your tree
unroot <- phy_tree(mag_phylo_sa)

# The circular layout tree
p <- ggtree(unroot, layout="fan", size=0.4, open.angle=0)

p
```

<img src="figures/ggtree to view unrooted tree-1.png" style="display: block; margin: auto;" />

#### Adding tip points with ggstar

```r
#add tip points with geom_star of ggstar
library(ggnewscale)

p1 <- p %<+% mag_full_df +
  geom_star(
    mapping=aes(fill=Phylum, starshape=fraction_class),
    starstroke=0.6
  ) +
  scale_fill_manual(
    values=c("#3f50a3", "#C0B956", "#e46009", "#CE8C1D","#965126","#553D32", "#04d2d5", "#80405e", "#b2e902", "#72c998", "#9FB6D9", "#B994E5", "#DF3383"),
    guide=guide_legend(keywidth = .5, keyheight = .5, 
                       order = 1, 
                       override.aes = list(starshape=15)),
                       na.translate=FALSE 
  ) +
  scale_starshape_manual(
    values=c(15,1,12),
    guide=guide_legend(keywidth = 0.5, keyheight = 0.5, 
                      labs(fill="Fraction Class"), order = 2),
    na.translate=FALSE
  ) +
  new_scale_fill() 

#p1 has ggstar that colors by bacterial phylym and has node tip points as fraction_class
p1
```

<img src="figures/tip points with ggstar-1.png" style="display: block; margin: auto;" />

```r
#notes about the clade info:
    #cyanobacteria: node 1-19, 348, 349(partial), 350(partial), 351-352, 353(partialbar), 354(partial), 355-365
    #Armatinomondota: node 21-22, maybe bar in node 367 and 368 need to adjust to see
    #Firmicutes node:23-27; node bar: 366 but also might be, 369 (thisisthebar) 370-372
    #Actinobacteriota node:27-63, 376-377(partialbar), 378-382, 383-387(partialbar), 388-391, 392(partial), 393, 394(partial), 395-398, 399-400(partial), 401, 402(partial), 403-409, 410-411(partial), 412, 413 (partial), 414-418, node bar: 374 and 375
    #Gemmatimonadota node:74-76, 422-423
    #Bacteroidota node: 77-155, 373, 421,424-501, node bar: 421 and 424 and 425
      #somany partials between 300-501
    #Bedellovibrionota and maybe _C: 339-344, 685-689 node bar 685 maybe node bar 684 or even 686
    #Bdellovibrionota or maybe _C: node bar 688
    #Verrucomicrobiota node: 156-191, 503-537,  node bar: 503
      #partial: 506, 509, 510 or 512, 513, 514-516, 525, stop checking track after 525 thres a lot (idk lazy to see)
    #Planctomyceotota node: 193-252, 538-597, node bar: 538
    #Proteobacteria node: 255-338, 601-683, node bar: 601
    #Acideobacteriota node: 253-254, node bar MAYBE: 599
    #Myxococcota_A node: 345-346, node bar: 690
    
    #extra: this wrapsaroundn like 90% of tree for some reason geom_cladelab(node = 373 , label = "wrapsaroundmostoftree")+
    #extra: node 347 wraps around entire tree
    #extra: node 419 wraps around like 75%ish of tree
    #extra: node 420 wraps around like 50% of bottom half of tree 
    #extra: this includes Verrucomicrogiota until Planctomyceotota geom_cladelab(node = 502 , label = "Verrucomicrobiota+ Planctomycetota")
    #extra: node 543 covers like 95% of Planctomyceotota but 538 is better geom_cladelab(node = 543 , label = "englishaccent")+
    #extra" node 598 and 600 covers Proteobacteria+Myxoccocta_A+ Bdellovironota+Bdellovibrionota_C +Acideobacteriota
    #extra: node 684 seems to cover Bdellovibrionota + Bdellovribrionota_C + Myxococcota_A but i ened to see more
```

#### Color tree with geom_tippoint instead of geom_star but ended up not using this for my own usage 

```r
# p1.0 <- p %<+% mag_full_df +
#   geom_tippoint(aes(color= fraction_class))+
#   scale_fill_manual(
#     values = c("#965953", "black", "#53977a"))
# 
# #p1 has ggstar that colors by bacterial phylym and has node tip points as fraction_class
# p1.0
```

#### Now using geom_hilight to highlight ALL bacterial phyla 

```r
#highlighting nodes
p2 <- p1 +
  geom_hilight(node = 599, alpha=0.4, fill="#3f50a3", extendto= 1.45, size=0.05) + #acido
  geom_hilight(node = 374, alpha=0.45, fill="#C0B956", extendto= 1.45, size=0.05) + #actino
  geom_hilight(node = 367, alpha=0.5, fill="#e46009", extendto= 1.45, size=0.05) + #arma
  geom_hilight(node = 424, alpha=0.4, fill="#CE8C1D", extendto= 1.45, size=0.05) + #bacterio
  geom_hilight(node = 685, alpha=0.5, fill="#965126", extendto= 1.45, size=0.3) + #bdello
  geom_hilight(node = 689, alpha=0.55, fill="#553D32", extendto= 1.45, size=0.3) + #bdello_c
  geom_hilight(node = 348, alpha=0.25, fill="#04d2d5", extendto= 1.45, size=0.05) + #cyano
  geom_hilight(node = 369, alpha=0.4, fill="#80405e", extendto= 1.45, size=0.05) + #firmicutes
  geom_hilight(node = 422, alpha=0.35, fill="#b2e902", extendto= 1.45, size=0.3) + #gemma
  geom_hilight(node = 690, alpha=0.55, fill="#72c998", extendto= 1.45) + #myxo
  geom_hilight(node = 538, alpha=0.4, fill="#9FB6D9", extendto= 1.45, size=0.05) + #planctomycetota
  geom_hilight(node = 601, alpha=0.4, fill="#B994E5", extendto= 1.45, size=0.05) + #proteo
  geom_hilight(node = 503, alpha=0.4, fill="#DF3383", extendto= 1.45, size=0.5)  #verru
  
p2
```

<img src="figures/geom_hilight-1.png" style="display: block; margin: auto;" />

#### Labeling clades with geom_cladelabel

```r
#labeling clades
#I only labeled the bigger phyla, the smaller ones i excluded for aesethetic purposes
p3 <- p2 +
  geom_cladelabel(node = 374, label= "Actinobacteriota", align= TRUE, offset.text=.1,barsize=NA, fontsize=5, angle="auto", colour = "#8a863e", hjust="center", horizontal=FALSE) +
  geom_cladelabel(node = 424, label= "Bacteriodota", align= TRUE, offset.text=.1,barsize=NA, fontsize=5, angle="auto", colour = "#a67117", hjust="center", horizontal=FALSE) + 
  geom_cladelabel(node = 348, label= "Cyanobacteria", align= TRUE, offset.text=.07,barsize=NA, fontsize=4, angle="auto", colour = "#03bdbf", hjust="center", horizontal=FALSE) +
  geom_cladelabel(node = 538, label= "Planctomycetota", align= TRUE, offset.text=.07,barsize=NA, fontsize=5, angle="auto", colour = "#8fa3c3", hjust="center", horizontal=FALSE) +
  geom_cladelabel(node = 601, label= "Proteobacteria", align= TRUE, offset.text=.07,barsize=NA, fontsize=5, angle="auto", colour = "#a685ce", hjust="center", horizontal=FALSE) + 
  geom_cladelabel(node = 503, label= "Verrucomicrobiota", align= TRUE, offset.text=.07,barsize=NA, fontsize=5, angle="auto", colour = "#DF3383", hjust="center", horizontal=FALSE) 

p3 #this has all the clade highlights and the larger clades are labelled
```

<img src="figures/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#heatmap of conjugation genes

```r
#heatmap, using geom_fruit to link geom_tile
#joining total_conj_gene_per_bin with mag_full using full_join to match everything tgt
  #be mindful that the count total_conj
testing_join <- full_join(congenee_all, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
#replacing NA in count genes
testing_join$count_conjugation[is.na(testing_join$count_conjugation)] <- 0

#first lets clean up our file because we keep getting errors that the columns are the same
new_total <- select(testing_join, -c(fraction_class.x,log2FoldChange,padj,fraction_class.y,Domain,Phylum,Class,Order,Family,Genus,Species,estimated_length.x))
```

```
## Error in `select()`:
## ! Can't select columns that don't exist.
## ✖ Column `log2FoldChange` doesn't exist.
```

```r
#convert new_total that is currently a tibble into dataframe
new_total_df <- as.data.frame(new_total)
```

```
## Error in eval(expr, envir, enclos): object 'new_total' not found
```

```r
#now lets make tree with heatmap ring with ggtreeExtra
gg_conj <- p2 + 
  geom_fruit(
    data=congenee_type,
    geom=geom_tile,
    mapping=aes(y=Seq_Name, x=conj_gene_type, alpha=conj_genes_per_bin, fill=conj_gene_type),
    color = "#525252",
    size = 0.35,
    offset = 0.05,
    pwidth=0.4
  ) +
  scale_alpha_continuous(
    range=c(0,5.5), 
    guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order= 5)
  ) +
  scale_fill_manual(
    values = c("#313167", "#aa5c3a", "olivedrab1", "#eba2b7"),
    guide=guide_legend(keywidth = 0.5, keyheight = 0.3, order = 4),
    na.translate= FALSE
  ) 
gg_conj #this is trying to look at a heatmap of the conjugation genes without the clade labels 
```

<img src="figures/ggtreeExtra conjugation genes heatmap-1.png" style="display: block; margin: auto;" />

#### geom_col to look at plasmid genes per bin 

```r
#we want to merge these numbers together and count them so that we can plot them as a bar graph!

testing_plasmid <- full_join(plasmid_sum_magfull_gene_count, mag_full, by = c("Seq_Name" = "Seq_Name"))
#replacing NA in count genes, make it numeric
testing_plasmid$count_n_genes[is.na(testing_plasmid$count_n_genes)] <- 0

#now lets try to get estimated genes per bin unsure if this is necessary 
testing_plasmid$genes_per_bin <- testing_plasmid$count_n_genes / testing_plasmid$estimated_length.y

#make it manipulatitable as data frame
testing_plasmid <- as.data.frame(testing_plasmid)

plasmid_sum_magfull_gene_count1 <- plasmid_sum_magfull_gene_count %>% 
  select(c("Seq_Name", "count_n_genes", "genes_per_bin", "fraction_class")) 

p4 <- p2 + 
  geom_fruit(
    data=plasmid_sum_magfull_gene_count,
    geom=geom_col,
    mapping=aes(x=genes_per_bin, y=Seq_Name, fill=fraction_class.y),
    pwidth=.5,
    position=position_stackx()
  ) +
  scale_fill_manual(
    values = c("#965953", "black", "#53977a")
  ) +
  geom_treescale(fontsize=2.9, linesize = .2, x= 2)
 
p4 #this is showing us the number of plasmid encoded genes that are divided by MAG
```

<img src="figures/geom_col plasmid genes per bin ggtreeExtra-1.png" style="display: block; margin: auto;" />

#### Phage Encoded Genes geom_col

```r
#using the phab_genes df from earlier


#bar graph of phage encoded genes. we will be using the dataframe count_phab_seq to see the number of plasmid encoded genes per bin that has the calculations for # phage genes per MAG
test_merge_phage <- full_join(count_phab_seq, counting_phage_genes, by=c("Seq_Name"="bin"))
```

```
## Error in eval(expr, envir, enclos): object 'count_phab_seq' not found
```

```r
#now turning the NAs into zeros
test_merge_phage$genes_per_bin[is.na(test_merge_phage$genes_per_bin)] <- 0
```

```
## Error: object 'test_merge_phage' not found
```

```r
#cleaning up our merged file
phab_genes1 <- phab_genes %>% 
  select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "CUBHE", "GC", "GCdiv", "ConsistencyHE", "CUB", "CPB", "nHE", "dCUB", "d", "LowerCI", "UpperCI", "FilteredSequences", "growth_class", "fraction_class.x", "Bin", "log2FoldChange","padj", "label", "angle", "Completeness", "fraction_class", "growth_class", "bin_length", "length_mb")) 
```

```
## Error in `select()`:
## ! Can't select columns that don't exist.
## ✖ Column `Domain` doesn't exist.
```

```r
phab_genes <- as.data.frame(phab_genes)
#now lets make the ring                         
phage_tree <- p2 +
  geom_fruit(
    data=phab_genes,
    geom=geom_col,
    mapping=aes(y=Seq_Name, x=genes_per_bin_phage, fill=fraction_class.y),
    orientation="y",
    pwidth=.38,
    position=position_stackx()
  ) +
  scale_fill_manual(
    values = c("#965953", "black", "#53977a")
  ) +
  geom_treescale(fontsize=2.9, linesize = .2, x=2)
 
phage_tree #this is showing us the number of plasmid encoded genes that are divided by MAG
```

<img src="figures/phage encoded genes heatmap-1.png" style="display: block; margin: auto;" />

```r
#now lets use ggsave to save this plot at a higher dpi and size
ggsave("phage_genes_phylotree.png", plot = phage_tree, width = 5, height = 4, dpi = 1000)
```

#### Pagels Lambda to describe trait evolution 

```r
#first we need to define which trait we want to test and give names to each value according to species
    #1 for example we could look at the trait of having conjugation genes 
trait <- testing_join
trait_conj <- testing_join$count_conjugation
names(trait_conj) <- testing_join$Seq_Name

unrooted <-phy_tree(mag_phylo_sa)
tree_label <- unrooted %>%
  as_tibble()%>%
  left_join(testing_join, by = c("label" = "Seq_Name"))%>%
  as.treedata()

#pagel's
tree_phylo <- tree_label %>% 
  as_tibble %>% as.phylo
phylosig(tree_phylo, trait_conj, method = "lambda", test=T)
```

```
## 
## Phylogenetic signal lambda : 0.992932 
## logL(lambda) : -512.203 
## LR(lambda=0) : 240.196 
## P-value (based on LR test) : 3.56436e-54
```

```r
# Phylogenetic signal lambda : 0.992932 
# logL(lambda) : -512.203 
# LR(lambda=0) : 240.196 
# P-value (based on LR test) : 3.56436e-54

#Bloomberg's K
phylosig(tree_phylo, trait_conj, method = "K", test=T, nsim=999)
```

```
## 
## Phylogenetic signal K : 0.0326016 
## P-value (based on 999 randomizations) : 0.001001
```

```r
# Phylogenetic signal K : 0.0326016 
# P-value (based on 999 randomizations) : 0.001001



#now lets see pagels and bloombergs but on plasmid encoded genes per bin

#using the testing_plasmid df we created earlier
trait_plasmid <- testing_plasmid$genes_per_bin
names(trait_plasmid) <- testing_plasmid$Seq_Name

unrooted1 <-phy_tree(mag_phylo_sa)
tree_label <- unrooted1 %>%
  as_tibble()%>%
  left_join(testing_plasmid, by = c("label" = "Seq_Name"))%>%
  as.treedata()

#pagels test
tree_phylo <- tree_label %>% 
  as_tibble %>% as.phylo
phylosig(tree_phylo, trait_plasmid, method = "lambda", test=T)
```

```
## 
## Phylogenetic signal lambda : 0.741024 
## logL(lambda) : 90.7589 
## LR(lambda=0) : 199.243 
## P-value (based on LR test) : 3.05552e-45
```

```r
# Phylogenetic signal lambda : 0.741024 
# logL(lambda) : 90.7589 
# LR(lambda=0) : 199.243 
# P-value (based on LR test) : 3.05552e-45 

#Bloomberg's K
phylosig(tree_phylo, trait_plasmid, method = "K", test=T, nsim=999)
```

```
## 
## Phylogenetic signal K : 0.00532172 
## P-value (based on 999 randomizations) : 0.002002
```

```r
# Phylogenetic signal K : 0.00532172 
# P-value (based on 999 randomizations) : 0.002002 


#now we can also see this for phage encoded genes
phage_presence_test <- tax_phass %>%
  group_by(Seq_Name, estimated_length) %>%
  summarise(
    count_virus_hits = sum(n_genes)
  )
```

```
## Error in `group_by()`:
## ! Must group by variables found in `.data`.
## ✖ Column `estimated_length` is not found.
```

```r
#now we can divide number of hits by estimated length
phage_presence_test$abundance <- phage_presence_test$count_virus_hits / phage_presence_test$estimated_length
```

```
## Error in eval(expr, envir, enclos): object 'phage_presence_test' not found
```

```r
testing_phage <- full_join(phage_presence_test, mag_full_df, by = c("Seq_Name" = "Seq_Name"))
```

```
## Error in eval(expr, envir, enclos): object 'phage_presence_test' not found
```

```r
#replacing NA in count genes
testing_phage$count_virus_hits[is.na(testing_phage$count_virus_hits)] <- 0
```

```
## Error: object 'testing_phage' not found
```

```r
testing_phage$abundance[is.na(testing_phage$abundance)] <- 0
```

```
## Error: object 'testing_phage' not found
```

```r
#using the testing_plasmid df we created earlier
trait_phage <- testing_phage$abundance
```

```
## Error in eval(expr, envir, enclos): object 'testing_phage' not found
```

```r
names(trait_phage) <- testing_phage$Seq_Name
```

```
## Error in eval(expr, envir, enclos): object 'testing_phage' not found
```

```r
unrooted2 <-phy_tree(mag_phylo_sa)
tree_label <- unrooted2 %>%
  as_tibble()%>%
  left_join(testing_phage, by = c("label" = "Seq_Name"))%>%
  as.treedata()
```

```
## Error: object 'testing_phage' not found
```

```r
tree_phylo_phage <- tree_label %>% 
  as_tibble %>% as.phylo
phylosig(tree_phylo_phage, trait_phage, method = "lambda", test=T)
```

```
## Error in eval(expr, envir, enclos): object 'trait_phage' not found
```

```r
# Phylogenetic signal lambda : 0.242194 
# logL(lambda) : -370.64 
# LR(lambda=0) : 6.77843 
# P-value (based on LR test) : 0.00922661 

#Bloomberg's K
phylosig(tree_phylo_phage, trait_phage, method = "K", test=T, nsim=999)
```

```
## Error in eval(expr, envir, enclos): object 'trait_phage' not found
```

```r
# Phylogenetic signal K : 0.00182091 
# P-value (based on 999 randomizations) : 0.927928 


#this is just from the data that we have not the WHOLE mag
trait_phage <- testing_phage$abundance
```

```
## Error in eval(expr, envir, enclos): object 'testing_phage' not found
```

```r
names(trait_phage) <- testing_phage$Seq_Name
```

```
## Error in eval(expr, envir, enclos): object 'testing_phage' not found
```

```r
testing_full_phage_merge <- left_join(phage_presence, mag_full_df, by = c("Seq_Name" = "Seq_Name"))

testing_full_phage_merge$Phylum_vir[is.na(testing_full_phage_merge$Phylum_vir)] <- 0

trait_phage_vir <- testing_full_phage_merge$Phylum_vir
names(trait_phage_vir) <- testing_full_phage_merge$Seq_Name

unrooted3 <-phy_tree(mag_phylo_sa)
tree_label <- unrooted3 %>%
  as_tibble()%>%
  left_join(testing_full_phage_merge, by = c("label" = "Seq_Name"))%>%
  as.treedata()
```

```
## Error in check_edgelist(x): Cannot find root. network is not a tree!
```

```r
tree_phylo_phage_vir <- tree_label %>% 
  as_tibble %>% as.phylo
phylosig(tree_phylo_phage, trait_phage_vir, method = "lambda", test=T)
```

```
## [1] "some species in tree are missing from x , dropping missing taxa from the tree"
```

```
## Error in invCl %*% y: requires numeric/complex matrix/vector arguments
```

```r
# Phylogenetic signal lambda : 0.242194 
# logL(lambda) : -370.64 
# LR(lambda=0) : 6.77843 
# P-value (based on LR test) : 0.00922661 

#Bloomberg's K
phylosig(tree_phylo_phage, trait_phage, method = "K", test=T, nsim=999)
```

```
## Error in eval(expr, envir, enclos): object 'trait_phage' not found
```

```r
# Phylogenetic signal K : 0.00182091 
# P-value (based on 999 randomizations) : 0.906907 
```

### Using bakta to annontate MAGs
So geNomad weirdly did not annontate antibiotic resistance genes so I wanna see if theres really nothing or if there are some hiding? It could be that the dataset that we are working with are short reads so maybe its having trouble annontating them because of that?


Anyways so today we are going to use bakta to annontate MAGs #thoughts and prayers

We will be running bakta through a conda environment!

At first I was planning to run through biohpc but their software wasnt updated and it was faster for me to work with a conda environment. I also had dependency issues once the biohpc updated the software so the conda environment was easier and is acutally working!

```bash
#so to run bakta we will need to download their database in order to annontate the DNA sequences we will make bakta folder inside our data folder and this is where we will run the analysis
cd /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/346_MAGs_bakta

#create conda environment
conda create --name bakta
conda activate bakta

#install with bioconda (takes a while but worth it pip was super quick but I also had dependency issues?)
conda install -c conda-forge -c bioconda bakta

#now lets download the database
  wget https://zenodo.org/record/10522951/files/db-light.tar.gz
  tar -xzf db-light.tar.gz
  rm db-light.tar.gz

#move to the correct directory that has the 346 MAG .fa files (I made a copy of the original files and put them into a folder called 346_MAGs_bakta)
cd Data/346_MAGs/346_MAGs_bakta
mkdir bakta_annontation

#my amrfinder db was weird and not workign so i forced it to update (omg because it was literally updated last week kms but yay for recent update!)
amrfinder_update --force_update --database bakta_database/amrfinderplus-db/
#running code with with bakta conda environment in a loop 
for file in *.fa
  do 
    tag=${file%.fa}  
    echo "Annotating $file with bakta for AMR genes"
    #--db is where the database is located
    bakta --skip-trna --skip-tmrna --skip-ncrna --skip-ncrna-region --skip-crispr --skip-sorf --skip-ori --meta --verbose --db bakta_database --prefix "$tag" --skip-trna --threads 40 --output bakta_annontation/"$tag"_bakta $file
    echo "Put output in data/MAG_bakta_annotation/'$tag'_bakta"
  done

echo "Prokka annotation finished!"

# verbose output writing results to results directory with  file prefix and  locus tag using an existing prodigal training file
```

### bakta AMR genes parsing
ok so bakta annontated all the CDSs so now we need to parse out just the AMR genes

This step will be annoying because it does not put the amr genes into a single column so it might have to be a manual thing if it is kill me now:)

```r
DAStool_91614_bin3 <- read.delim("91614_DAStool_bin3.tsv")
```

```
## Error in file(file, "rt"): cannot open the connection
```

### bakta: so its bizarre to me that it does not find any AMR genes even though I CLEARLY see them being annontated like i literally see something with vancomycin resistance and a bunch of beta lactamases so I think maybe its just not liking the fact I have metagenomic sequences? so I will try another program just to be sure. I did try downloading the nt seqs and uploading one with a clear arg that I see and something that geNomad picked up and I ran it through resfinder but it didnt annontate it with any resistance at all which is nuts. and i see that there are no pseudo genes in those files from the log produced by bakta so i do not know what the heck is going on. I also tried staramr but i realized after teh fact when it wasnt working, its only meant for whole genome assemblies with a specific organism in mind, not metagenomic samples/ assembled genomes. 
### thsi paper (https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0648-z) used a program called ABRICATE to annontate their ARGs and virulence genes. their paper also had MAGs  so maybe it will work?


```bash
#install abricate
conda create --name abricate

conda install -c conda-forge -c bioconda -c defaults abricate
abricate --check
abricate --list

#make sure you are in the current working directory with 346 MAG *.fa files
abricate -db ncbi --minid 70 --mincov 70 *.fa > abricate_results.txt

#i was able to get one gene annontated with teh 70% cutoff and this was also validated by NCBI AMRfinder!! weirdly not by geNomad but this is a win
```

### just for fun adn to see i found another program that might be good for annontating metagenomic ARGS its called fargene. here is a snippet of whats on their github!

"fARGene (Fragmented Antibiotic Resistance Gene iENntifiEr ) is a tool that takes either fragmented metagenomic data or longer sequences as input and predicts and delivers full-length antiobiotic resistance genes as output"

```bash
#create condda environment
conda create --name fargene
conda activate fargene
#install from bioconda
conda install -c conda-forge -c bioconda fargene

#cd to where files are
cd Data/fargene

#example of how to run
fargene -i path/to/fastafile(s)/*.fasta --no-quality-filtering --no-orf-prediction --hmm-model class_a -o output_dir

fargene -i path/to/fastafile(s)/*.fasta --hmm-model class_a -o output_dir

#fargene is an older program and it has been been updated since 2022 but it still works fine overall. like half the genes that you can search for work automatically. there are some cases like the aminoglycoside genes or anything that throws and error because those genes and their HMM profiles are not added in the conda environment for fargene. to get around that you must install the conda package AND clone the github repo into your folder and provide a path for the gene you are interested in. it is more manual but it isnt bad whatsoever. if anyone were to repeat this you could follow this code and be fine. be aware that the fargene hmm path only wants it in an absolute path which is a little silly especially for reproducibility but other than that it works! because it si more manual you will need to manually enter the --score and --meta score these values were found on the cloned git repo in a python folder. for more clarification there is an issue on the fargene page. 

#beta lactam

#class a
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --no-quality-filtering --no-orf-predict --hmm-model class_a -o class_a_out -p 30
#class b_1_2
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model class_b_1_2 -o class_b_out -p 30 
#class_b_3
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model class_b_3 -o class_b_3_out -p 30 
#class_c
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model class_c -o class_c_out -p 30 
#class_d_1
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model class_d_1 -o class_d1_out -p 30 
#class_d_2
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model class_d_2 -o class_d2_out -p 30 

#b1
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/B1.hmm -o B1 --score 135.8 --meta-score 0.2424 -p 30

#tetracycline

#tet efflux
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model tet_efflux -o tet_efflux_out -p 30 
#tet rpg
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model tet_rpg -o tet_rpg -p 30 
#tet enzyme
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa  --hmm-model tet_enzyme -o tet_enzyme -p 30 

#mph
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/mph.hmm -o mph --score 100 --meta-score 0.3030 -p 30

#qnr
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model qnr -o qnr_out -p 30 

#aminoglycoside 
#a
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/aminoglycoside_model_a.hmm -o aminogly_a --score 100 --meta-score .2727 -p 30
#c
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/aminoglycoside_model_c.hmm -o aminogly_c --score 100 --meta-score 0.2424 -p 30
#f
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/aminoglycoside_model_f.hmm -o aminogly_f1 --score 100 --meta-score .3636 -p 30
#h
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/aminoglycoside_model_h.hmm -o aminogly_h --score 100 --meta-score .1818 -p 30
#i
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/aminoglycoside_model_i.hmm -o aminogly_i --score 100 --meta-score 0.2121 -p 30

#erm
#type a
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/erm_typeA.hmm -o erm_typeA --score 200 --meta-score 0.4545 -p 30
#type f
fargene -i /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/*.fa --hmm-model /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/fargene/fargene/fargene_analysis/models/erm_typeF.hmm -o erm_typeF --score 175 --meta-score 0.4242 -p 30
```

## Visualizing ARGs

so we will first get our ARGs that have been annontated by geNomad in the amr_gene_file column. the annoying thing is that its in database format so idk what it is easily but thats fine we can fix!

we will be using casewhen and creating one column for ncbi gene name id, then we will be translating the ncbi name id into a readable resistance gene, then another column for antibiotic class, and maybe another column for resistance mechanism?


```r
#for the ncbi name im omitting the part that says -NCBI so really for example the full name is Fos_GSH-NCBIFAM but im just doing Fos_GSH
df_plasmidsum_clean_up <- df_plasmidgenes %>% 
  mutate(ARG_clean = case_when(annotation_amr == "NF000088" ~ "viomycin_Vph",
                               annotation_amr == "NF033135" ~ "cmx_cmrA",
                               annotation_amr == "NF012178" ~ "tet_MFS_V",
                               annotation_amr == "NF000496" ~ "Fos_GSH",
                               annotation_amr == "NF033117" ~ "vanR_ACDEGLN",
                               annotation_amr == "NF033476" ~ "tet_destruct",
                               annotation_amr == "NF033134" ~ "cmlA_floR",
                               annotation_amr == "NF033156" ~ "ble_BLMA_gen",
                               annotation_amr == "NF033693" ~ "AAC_6p_Ie_fam",
                               annotation_amr == "NF033105" ~ "bla_subclass_B3",
                               annotation_amr == "NF000165" ~ "AAC_6p_Entco",
                               annotation_amr == "NF000272" ~ "OqxA_adapt_sub",
                               annotation_amr == "NF012171" ~ "APH_6",
                               annotation_amr == "NF012174" ~ "tet_MFS_A_B_C_D",
                               annotation_amr == "NF000402" ~ "vanR-B",
                               annotation_amr == "NF033088" ~ "bla_subclass_B1",
                               annotation_amr == "NF000005" ~ "ble_BLMT",
                               is.na(annotation_amr) ~ "NA"))

# now we will create a new column detailing the resistance mechansim of the proposed gene. The mechanisms are based off of the card database (https://card.mcmaster.ca/ontology/39137)
df_plasmidsum_clean_up <- df_plasmidsum_clean_up %>% 
  mutate(resistance_mech = case_when(ARG_clean == "viomycin_Vph" ~ "Antibiotic Inactivation",
                               ARG_clean == "cmx_cmrA" ~ "Antibiotic Efflux",
                               ARG_clean == "tet_MFS_V" ~ "Antibiotic Efflux",
                               ARG_clean == "Fos_GSH" ~ "Antibiotic Inactivation",
                               ARG_clean == "vanR_ACDEGLN" ~ "Antibiotic Target Alteration",
                               ARG_clean == "tet_destruct" ~ "Antibiotic Inactivation",
                               ARG_clean == "cmlA_floR" ~ "Antibiotic Efflux",
                               ARG_clean == "ble_BLMA_gen" ~ "Antibiotic Inactivation",
                               ARG_clean == "AAC_6p_Ie_fam" ~ "Antibiotic Inactivation",
                               ARG_clean == "bla_subclass_B3" ~ "Antibiotic Inactivation",
                               ARG_clean == "cmx_cmrA" ~ "Antibiotic Efflux",
                               ARG_clean == "AAC_6p_Entco" ~ "Antibiotic Inactivation",
                               ARG_clean == "OqxA_adapt_sub" ~ "Antibiotic Efflux",
                               ARG_clean == "APH_6" ~ "Antibiotic Inactivation",
                               ARG_clean == "tet_MFS_A_B_C_D" ~ "Antibiotic Efflux",
                               ARG_clean == "vanR-B" ~ "Antibiotic Target Alteration",
                               ARG_clean == "bla_subclass_B1" ~ "Antibiotic Inactivation",
                               ARG_clean == "ble_BLMT" ~ "Antibiotic Inactivation",
                               is.na(ARG_clean) ~ "NA"))

# now id like to have a new column that will have a better explanation of each of the antibiotics and the types of antibiotics. so if i have a bla gene it will be classified under "beta lactam"
df_plasmidsum_clean_up <- df_plasmidsum_clean_up %>% 
  mutate(antibiotic_name = case_when(ARG_clean == "viomycin_Vph" ~ "Viomycin",
                               ARG_clean == "cmx_cmrA" ~ "Chloramphenicol",
                               ARG_clean == "tet_MFS_V" ~ "Tetracycline",
                               ARG_clean == "Fos_GSH" ~ "Fosfomycin",
                               ARG_clean == "vanR_ACDEGLN" ~ "Vancomycin",
                               ARG_clean == "tet_destruct" ~ "Tetracycline",
                               ARG_clean == "cmlA_floR" ~ "Chloramphenicol",
                               ARG_clean == "ble_BLMA_gen" ~ "Bleomycin",
                               ARG_clean == "AAC_6p_Ie_fam" ~ "Aminoglycoside",
                               ARG_clean == "bla_subclass_B3" ~ "Beta Lactam",
                               ARG_clean == "cmx_cmrA" ~ "Chloramphenicol; Bleomycin",
                               ARG_clean == "AAC_6p_Entco" ~ "Aminoglycoside",
                               ARG_clean == "OqxA_adapt_sub" ~ "Fluoroquinolone",
                               ARG_clean == "APH_6" ~ "Aminoglycoside",
                               ARG_clean == "tet_MFS_A_B_C_D" ~ "Tetracycline",
                               ARG_clean == "vanR-B" ~ "Vancomycin",
                               ARG_clean == "bla_subclass_B1" ~ "Beta Lactam",
                               ARG_clean == "ble_BLMT" ~ "Bleomycin",
                               is.na(ARG_clean) ~ "NA"))

#now to get to the visualization part lets merge this dataset with the magfull to see where the ARGs fall in each of the bacterial lifestyles
plasmid_geno_args_magfull <- left_join(df_plasmidsum_clean_up, mag_full_df, by = c("Seq_Name" = "Seq_Name"))


# ok now we have merged the documents and we can do the fun part of visualizing the genes, their fraction class, type of antibiotic 

#to make a bar graph we are going to need to count the number of genes
#filter out nas
plasmid_geno_args_magfull_no_na <- plasmid_geno_args_magfull[!is.na(plasmid_geno_args_magfull$annotation_amr), ]

count_args_geno <- plasmid_geno_args_magfull_no_na %>%
  group_by(fraction_class, antibiotic_name ) %>%
  summarise(count_args = n())

args_geNo <- ggplot(count_args_geno, aes(x = fraction_class, y = count_args, fill = antibiotic_name))+
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Bacterial Lifestyle")+
  ylab("Antibiotic Resistance Genes")+
  scale_fill_wa_d(palette = "rainier", reverse = TRUE) +
  labs(fill = "args") +
  theme_classic()
args_geNo
```

<img src="figures/geNomad-ARGs-clean-up-1.png" style="display: block; margin: auto;" />

```r
#now lets use ggsave to save this plot at a higher dpi and size
ggsave("geNomad-ARGs-clean-up.png", plot = args_geNo, width = 5, height = 4, dpi = 1000)
```

#### fargene analysis

```r
#class a beta lactam
  # when I look at the predicted orfs nt and amino acid seqs I see that there are 3 predicted orfs when i put the nt seqs into blastx (translated nucleotides) i see that these truly are class a beta lactams. the only thing is technically i see that there ar e10 genes that are predicted but no open reading frame so im just gonna go with the 3 predicted orfs
  # when i blast I am getting pretty confident scores that this is class a beta lactam like above 80% query and percent identity

  # need to add 91616_bin_48_k99_1416201_1,free-coassembly-bin374, particle-co-bin374

#class b_1_2 beta lactam
  #2 predicted genes 0 orfs
  # need to add 91623_bin7, 91627_bin136 (had query + identity above 85%)

#class b_3 beta lactam
  # 12 predicted genes 9 predicted ORFs

## ok silly of me I accidentally deleted my analysis with the manual models I will add them in here. I have detailed notes of what had genes and what didnt so I will only 


# fargene analysis
```

#### log2fold of how mags were categorized into lifestyles
for my poster the background of how we categorized MAGs into lifestyles is needed for context. we will make a box plot or bar graph of the log2fold change where it will have shapes that denote lifestyle and will be separated by bacterial phyla levels

```r
#lets reorder our samples to have it be read like free living, particle attached, and then generalist
mag_full_df$fraction_class <- factor(mag_full_df$fraction_class, levels = c("Free", "Generalist", "Particle"))

#lets establish fraction class colors 
fraction_class_colors <- c(
  "Free" = "#C3774D", 
  "Generalist" = "#932252",
  "Particle" = "#8B9322")

#now lets visualize
log2fold_mags <- ggplot(mag_full_df, aes(x = log2FoldChange, y = Phylum, fill = fraction_class, shape = fraction_class)) +
  geom_col(position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.4), show.legend = TRUE,) +
  scale_fill_manual(values = fraction_class_colors) +
  scale_shape_manual(values = c("Free" = 24, "Generalist" = 21, "Particle" = 25)) +
  guides(shape = guide_legend(override.aes = list(size = 3.3))) +
  theme_classic() +
  theme(legend.position = "bottom") 
log2fold_mags
```

<img src="figures/log2fold-MAGs-1.png" style="display: block; margin: auto;" />

```r
ggsave("log2foldmag.png", plot = log2fold_mags, width = 5, height = 4, dpi = 1000)
```

