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


## this is me trying to understand the files so its a little chaotic sorry ##

#lets start with looking at proviruses
#set working directory
getwd()
setwd("/local/workdir/sna49/phage_personal_SA/genomad_practice/genomad_output/GCF_009025895.1.fna_find_proviruses/")

#some files are in .tsv formats so we are installing readr to help us analyze the results
install.packages("readr")
library(readr)

#lookng at provirus taxonomy
provirus_taxonomy <- read.table(file = "GCF_009025895.1.fna_provirus_taxonomy.tsv", sep ="\t", header = TRUE)
View(provirus_taxonomy)

#now lets look at the summary files

#1. looking at plasmid genes
setwd("/local/workdir/sna49/phage_personal_SA/genomad_practice/genomad_output/GCF_009025895.1.fna_summary")
plasmid_genes <- read.table(file = "GCF_009025895.1.fna_plasmid_genes.tsv", sep ="\t", header = TRUE)
View(plasmid_genes)

#2. looking at plasmid summary 
setwd("/local/workdir/sna49/phage_personal_SA/genomad_practice/genomad_output/GCF_009025895.1.fna_plasmid_summary")
plasmid_summary <- read.table(file = "GCF_009025895.1.fna_plasmid_summary.tsv", sep ="\t", header = TRUE)
View(plasmid_summary)

#3 looking at virus genes
setwd("/local/workdir/sna49/phage_personal_SA/genomad_practice/genomad_output/GCF_009025895.1.fna_summary")
virus_genes <- read.table(file = "GCF_009025895.1.fna_virus_genes.tsv", sep ="\t", header = TRUE)
View(virus_genes)

#4 looking at virus summary
setwd("/local/workdir/sna49/phage_personal_SA/genomad_practice/genomad_output/GCF_009025895.1.fna_summary")
virus_summary <- read.table(file = "GCF_009025895.1.fna_virus_summary.tsv", sep ="\t", header = TRUE)
View(virus_summary)
