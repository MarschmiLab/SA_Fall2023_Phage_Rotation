#this shell script will take geNomad output files from genomad_out_346 (my genomad output directory holding all files) into their respective summary files
#in summary it will take separate summary files and combine them inot a big summary output file

#set up environment to run genomad
#source /programs/miniconda3/bin/activate genomad
#genomad download-database .

pwd #make sure that you are in the correct path /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/genomad_output_346

#our results with the summary .tsv files are embedded in the output directory so we must get the .tsv files out of the directory. we do not want the log files.
#FIRST THINGS FIRST we want to get the files WITHIN the summary directories of each output.

#but we need to actually start by naming the actual directory where we want the results sent to
sum_directory="summary_tsv_files"

#now we will create the output directory
mkdir -p "sum_directory" #the -p ensures that there will be no error if there is an existing directory otherwise it will make the directory as needed

for d in *.fa_summary/ ; do #iterating thorugh directories ending in *.fa_summary/
  echo "$d" #printing each directory name with this ending
  find "$d" -type f -name "*.tsv" -exec cp {} "sum_directory/" \;
done

#now all files will be placed in ~/Ouput/genomad_output_346/sum_directory
#within that directory you will notice that all files have an ending with .tsv, but we will concatenate all files to end in their respective .txt file

cd sum_directory #make sure you are in sum_directory

#\\\\\\\wait maybe ignore this whole section kms////////
#1. lets concatenate ALL the plasmid_genes files
cat *.fa_plasmid_genes.tsv > sorted_plasmid_genes.txt
#2. lets concatenate ALL the plasmid_summary files
cat *.fa_plasmid_summary.tsv > sorted_plasmid_summary.txt
#3. lets concatenate ALL the virus_genes files
cat *.fa_virus_genes.tsv > sorted_virus_genes.txt
#4. lets concatenate ALL the plasmid_summary files
cat *.fa_virus_summary.tsv > sorted_virus_summary.txt
#\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////

#so now we will combine the output concatenated summary files into a long list and we will add a column that will specify the filename for each row

#so now we are using the awk command that processes certain patterns and actions
  #goes like: awk 'pattern { action }' input-file
  #i saw this in azealia's file and it seems useful and hopefully it works!


#\\\\\\\\\\\\ using awk with plasmid genes to create long list of plasmid genes with all the files //////////
#so you could make the .tsv files into a long list of .txt but I will be making them into .tsv files so I can analyze them in rstudio later (my earlier script already as .tsv files able to be read and analyzed)
#1. creating long list of all plasmid genes
awk 'NR==1{print $0" filename" }' *.fa_plasmid_genes.tsv > plasmid_genes.tsv
  #NR==1: checking to see if the current record number is 1, so it is processing first line aka the header of each file
  #print $0: filename: if the condition is met it will print the entire line ('$0') followed by space and string called "filename"
awk 'FNR>1{ sub(/\/[^\/]+$/, "", FILENAME); print FILENAME, $0 }' *.fa_plasmid_genes.tsv >> plasmid_genes.tsv
  #FNR>1: checks to see if current record number within the current file is greater than 1, skips the first line or the header of each input file
  #{ sub(/\/[^\/]+$/, "", FILENAME); print FILENAME, $0 }' : for reocrds that meet these conditions, it will modify the "filename" and remove teh directory path and then print modify 'filename' followed by entire line ('$0')
  
#2. creating long list of all plasmid summary
awk 'NR==1{print $0" filename" }' *.fa_plasmid_summary.tsv > plasmid_summary.tsv
awk 'FNR>1{ sub(/\/[^\/]+$/, "", FILENAME); print FILENAME, $0 }' *.fa_plasmid_summary.tsv >> plasmid_summary.tsv

#3. creating long list of all virus genes
awk 'NR==1{print $0" filename" }' *.fa_virus_genes.tsv > virus_genes.txt
awk 'FNR>1{ sub(/\/[^\/]+$/, "", FILENAME); print FILENAME, $0 }' *.fa_virus_genes.tsv >> virus_genes.tsv
 
#4. creating long list of all virus summary
awk 'NR==1{print $0" filename" }' *.fa_virus_summary.tsv > virus_summary.txt
awk 'FNR>1{ sub(/\/[^\/]+$/, "", FILENAME); print FILENAME, $0 }' *.fa_virus_summary.tsv >> virus_summary.tsv

#at this point we have successfully created our long lists of data of plasmid and phage information and we can now begin our Rstudio analysis



