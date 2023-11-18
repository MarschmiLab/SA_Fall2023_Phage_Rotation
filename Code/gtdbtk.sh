#downloading conda to computer
#I am writing how I was able to download conda and initialize it

#cd to workdir/"user" directory:
cd /local/workdir/sna49

#I used the quick command line to install
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

#after installing initailize your newly installed Miniconda. The following is the command for bash
~/miniconda3/bin/conda init bash

#then I was told to restart my shell for changes to take effect

#once installed, every time you need to use conda, run this command 
source $HOME/miniconda3/bin/activate

#///////////////////////////////////////////////////////////////////////////
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

#checking which version just in case
gtdbtk check_install

#gtdbtk does not like .fa.gz files just .fa so lets copy over the files 
cp -r /local/workdir/lab_data/2014-2015-metaG-muskegon/binning/dereplicated_bins /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk

#created new directory called gtdbtk 
/local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/346_MAGs
  #I recopied the original 346 MAGs but jsut the .fa version
  
#then go to gtdbtk_files and make symbolic link 
cd /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/gtdbtk_files/346_MAGs
#ln -s /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/gtdbtk_files/*.fa.gz
#ln -s /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/gtdbtk_files/346_MAGs

#update: I was not able to make a symbolic link and have it work for me so I just copied the files over

#now that we see that we have everything installed (yay!)
#we will now begin to do our analysis with the command called classify_wf 
  #this workflow consists fo 4 steps:
    #ani_screen: compares our genomes against a MASH database composed of all GTDB representative genoems then verifies the best mash hits using FASTANI which is a method used to estimate average nucleotide identity using alignment free approximate sequence mapping. it is good for finished and draft genomes
    #identify: uses prodigal and uses HMM mdoels and HMMER package to identify the 120 bacterial and 53 archaeal marker genes used for phylogenetic inference. MSA obtained by aligning amrker gene to respective HMM model
    #align: concatenates the aligned marker genes and filters the concatenated MSP to approx 5k amino acids
    #classify: uses pplacer to find the ML placement of each genome in teh gtdbtk reference tree. GTDBTK classifies each genome based on its placment in refernece tree, relative evolutionary divergence, and/or ANI to reference geomes
    
    #workflow command: gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
      #genomes must be in FASTA format but .gz is acceptable which is awesome bc thats what we got
      
      #### DISCLAIMER DO NOT USE THE CODE BELOW USE THE ONE THAT HAS --skip_ani_screen
# gtdbtk classify_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/346_MAGs --out_dir classify_gtdbtk_output2 --mash_db 2.3.2 -x .fa --cpus 5
  #gtdbtk = command
  #classify_wf = pipeline used
  ##--genome_dir = pointing to the genome directory or wherever our files are at
  #out_dir = pointing to the output directory
  #classify_gtdbtk_output = this is the output file that is being created
  #--mash_db 2.3.2 = this version of GTDBTK requires that you input the mash_db version (2.3.2 is the one we downloaded). This Using the --mash_db option will indicate to GTDB-Tk the path of the sketched Mash database require for ANI screening.
    #If no database are available ( i.e. this is the first time running classify ), the --mash_db option will sketch a new Mash database that can be used for subsequent calls.
  #-x .fa.gz = -x means extension files. the files we have are .fa.gz
  #--cpus 100 = these are the cpus that are being used by the server to calculate

#general notes about downloading
  #so far it is using Mash version 2.3
  #creating mash sketch file: classify_gtdbtk_output/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
  #completed 346 genomes in 1.00 seconds (347.67 genome)
  
    # [2023-10-21 00:04:17] INFO: GTDB-Tk v2.3.2
    # [2023-10-21 00:04:17] INFO: gtdbtk classify_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/gtdbtk_files --out_dir classify_gtdbtk_output --mash_db 2.3.2 -x .fa.gz --cpus 100
    # [2023-10-21 00:04:17] INFO: Using GTDB-Tk reference data version r214: /home/sna49/miniconda3/envs/gtdbtk-2.3.2/share/gtdbtk-2.3.2/db
    # [2023-10-21 00:04:18] INFO: Loading reference genomes.
    # [2023-10-21 00:04:19] INFO: Using Mash version 2.3
    # [2023-10-21 00:04:19] INFO: Creating Mash sketch file: classify_gtdbtk_output/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    # [2023-10-21 00:04:21] INFO: Completed 346 genomes in 1.31 seconds (263.36 genomes/second).
    # [2023-10-21 00:04:21] INFO: Creating Mash sketch file: 2.3.2.msh                  
    # [2023-10-21 00:08:02] INFO: Completed 85,205 genomes in 3.69 minutes (23,111.61 genomes/minute)

#it said that some of my files (specificslly the DASTool files were not gzipped correctly so I am getting the original files just in case the ones I had were corrupted somehow?)
  #cp -r /local/workdir/lab_data/2014-2015-metaG-muskegon/binning/dereplicated_bins/346_MAGs /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk
  
#so far GTDBTK gives us the msa alignment when we use the code above, this can be used for an unrooted tree which is totally fine


#now lets try classify_wf but with the --skip_ani_screen flag. This was recommended by the developers to get all 346 genomes characterized. The documentation does not reflect thtat BUT I am really happy it might work
gtdbtk classify_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/dereplicated_bins/346_MAGs --out_dir classify_gtdbtk_output_skip --skip_ani_screen -x .fa --cpus 5

#but we can also try the de novo wf to give us a rooted tree and we would use planctomycetes as our outgroup
gtdbtk de_novo_wf --genome_dir /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/gtdbtk/346_MAGs --outgroup_taxon p__Chloroflexota --bacteria  --out_dir de_novo_wf_test1 -x .fa --cpus 100




