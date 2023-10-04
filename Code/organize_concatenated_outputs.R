#use geNomad from bioHPC server

#before we begin we need to start the Rstudio server
#/programs/rstudio_server/mv_dir #this keeps the very big session data files under /workdir
#/programs/rstudio_server/rstudio_start #this is to start the rstudio server

#set up environment
#source /programs/miniconda3/bin/activate genomad

#run genomad
#genomad download-database .

#command to execute geNomad on shell
#genomad end-to-end --cleanup --splits 8 GCF_009025895.1.fna.gz genomad_output genomad_db

#once it runs there will be an output of a bunch of files 
#read .tsv files

#on rstudio



#lets start with looking at proviruses
#set working directory
getwd()
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")

#some files are in .tsv formats so we are installing readr to help us analyze the results
install.packages("readr")
library(readr)

rm(virus_summary)

#now lets look at the summary files

#1. looking at plasmid genes
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
plasmid_genes <- read.table(file = "plasmid_genes.tsv", sep ="\t", header = TRUE)
View(plasmid_genes)

#2. looking at plasmid summary 
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
plasmid_summary <- read.table(file = "plasmid_summary.tsv", sep ="\t", header = TRUE)
View(plasmid_summary)

#3 looking at virus genes
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
virus_genes <- read.table(file = "virus_genes.tsv", sep ="\t", header = TRUE)
View(virus_genes)

#4 looking at virus summary
setwd("/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/4_concatenated_files/")
virus_summary <- read.table(file = "virus_summary.tsv", sep ="\t", header = TRUE)
View(virus_summary)

#\\\\\\\\\\\\\\\\\\\\\\install or load library of required packages/////////////////////////
install.packages("tidyverse")
library(tidyverse)

install.packages("stringr")
library(stringr)

install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

install.packages("dplyr")
library(dplyr)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////
#1. plasmid genes
#need to create dataframe
df_plasmidgenes <- data.frame(plasmid_genes)

#split the first column called gene into multiple columns
#because there is so much going on with the column name we are going to separate each word into its own column and then mutate it together (dont be alarmed)
df_plasmidgenes <- df_plasmidgenes %>% separate(gene, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'plasmid', 'genes', 'tsv', 'Contig', 'bp_length', 'Contig_Number'))

#view to make sure all the outputs are looking good
View(df_plasmidgenes)

#now we are going to mutate the combine the columns together
df_plasmidgenes$Seq_Name <- paste(df_plasmidgenes$Seq_Name, df_plasmidgenes$DASTool, sep = "_")
df_plasmidgenes$fa_info <- paste(df_plasmidgenes$fa_info, df_plasmidgenes$plasmid, df_plasmidgenes$genes, df_plasmidgenes$tsv, sep = "_")
df_plasmidgenes$Contig <- paste(df_plasmidgenes$Contig, df_plasmidgenes$bp_length, sep = "_")

#view again to make sure it looks good
View(df_plasmidgenes)

#then we need to get rid of the extra columns that we do not need
df_plasmidgenes <- subset(df_plasmidgenes, select = -c(DASTool, plasmid, genes, tsv, bp_length)) #invalid argument need to fix

#view one more time to make sure it looks good
View(df_plasmidgenes)
#YAY IT LOOKS GREAT
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////
#2. plasmid summary
#need to create dataframe
df_plasmidsum <- data.frame(plasmid_summary)

#split the first column called gene into multiple columns
#because there is so much going on with the column name we are going to separate each word into its own column and then mutate it together (dont be alarmed)
df_plasmidsum <- df_plasmidsum %>% separate(seq_name, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'plasmid', 'genes', 'tsv', 'Contig', 'bp_length'))

#view to make sure all the outputs are looking good
View(df_plasmidsum)

#now we are going to mutate the combine the columns together
df_plasmidsum$Seq_Name <- paste(df_plasmidsum$Seq_Name, df_plasmidsum$DASTool, sep = "_")
df_plasmidsum$fa_info <- paste(df_plasmidsum$fa_info, df_plasmidsum$plasmid, df_plasmidsum$genes, df_plasmidsum$tsv, sep = "_")
df_plasmidsum$Contig <- paste(df_plasmidsum$Contig, df_plasmidsum$bp_length, sep = "_")

#view again to make sure it looks good
View(df_plasmidsum)

#then we need to get rid of the extra columns that we do not need
df_plasmidsum <- subset(df_plasmidsum, select = -c(DASTool, plasmid, genes, tsv, bp_length)) 

#view one more time to make sure it looks good
View(df_plasmidsum)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////
#3. virus genes
#need to create dataframe
df_virusgenes <- data.frame(virus_genes)

#split the first column called gene into multiple columns
#because there is so much going on with the column name we are going to separate each word into its own column and then mutate it together (dont be alarmed)
df_virusgenes <- df_virusgenes %>% separate(gene, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'virus', 'genes', 'tsv', 'Contig', 'bp_length', 'Contig_Number'))

#view to make sure all the outputs are looking good
View(df_virusgenes)

#now we are going to mutate the combine the columns together
df_virusgenes$Seq_Name <- paste(df_virusgenes$Seq_Name, df_virusgenes$DASTool, sep = "_")
df_virusgenes$fa_info <- paste(df_virusgenes$fa_info, df_virusgenes$virus, df_virusgenes$genes, df_virusgenes$tsv, sep = "_")
df_virusgenes$Contig <- paste(df_virusgenes$Contig, df_virusgenes$bp_length, sep = "_")

#view again to make sure it looks good
View(df_virusgenes)

#then we need to get rid of the extra columns that we do not need
df_virusgenes <- subset(df_virusgenes, select = -c(DASTool, virus, genes, tsv, bp_length)) #invalid argument need to fix

#view one more time to make sure it looks good
View(df_virusgenes)
#YAY IT LOOKS GREAT
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////
#4. virus summary
#need to create dataframe
df_virussum <- data.frame(virus_summary)

#split the first column called gene into multiple columns
#because there is so much going on with the column name we are going to separate each word into its own column and then mutate it together (dont be alarmed)
df_virussum <- df_virussum %>% separate(seq_name, c('Seq_Name', 'DASTool', 'Bin', 'fa_info', 'virus', 'genes', 'tsv', 'Contig', 'bp_length'))

#view to make sure all the outputs are looking good
View(df_virussum)

#now we are going to mutate the combine the columns together
df_virussum$Seq_Name <- paste(df_virussum$Seq_Name, df_virussum$DASTool, sep = "_")
df_virussum$fa_info <- paste(df_virussum$fa_info, df_virussum$virus, df_virussum$genes, df_virussum$tsv, sep = "_")
df_virussum$Contig <- paste(df_virussum$Contig, df_virussum$bp_length, sep = "_")

#view again to make sure it looks good
View(df_virussum)

#then we need to get rid of the extra columns that we do not need
df_virussum <- subset(df_virussum, select = -c(DASTool, virus, genes, tsv, bp_length)) 

#view one more time to make sure it looks good
View(df_virussum)

