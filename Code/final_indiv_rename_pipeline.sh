#starting with community analysis with geNomad with the final_indiv_assemblies
cp -r /local/workdir/lab_data/2014-2015-metaG-muskegon/assemblies/final_indiv_assemblies/ /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/final_indiv_assemblies/
#the -r means it will recursively copy the direcotry and its contents, if you do not include it shell gets mad and gives u error

#using the rename_geNomad.sh script to compress fasta files
#compressing all fasta files to run with geNomad

#all the files from the metagenomic data end in .fa and we are writing a script that will add .gz at the end of all the files

#make sure that you are in the correct directory, we should be in the /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/final_indiv_assemblies
for file in *.fa
do
  mv "$file" "${file%.fa}.fa.gz"
#gzip is needed to compress all the files
  gzip "${file%.fa}.fa.gz"
  echo "renamed files to .gz"
done

echo "all done!"

#/////////////////////////////////////////
#now we will run geNomad

#set up environment to run genomad
source /programs/miniconda3/bin/activate genomad
genomad download-database .

pwd #make sure that you are in the correct path /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/final_indiv_assemblies

#now time to run the loop to read the files and annonate them!
for file in *.fa.gz; 
  do 
  # Which file am I on? 
  echo $file 
  #genomad end-to-end = command that executes the pipeline of starting from nt FASTA file, then annonates, then finds proviruses, which then does marker classification to see if it is a chromosome, plasmid, or virus based on marker content, then does neural network (nn) classification but does not rely on marker infomration but processes it using IGLOO neural network, then does aggregated-classification to combine marker adn nn classifications, then does score calibration 
  #--cleanup = used to force geNomad to delete intermediate files that were generated to save space
  #--splits 32 = to prevent the execution form seaching in database of protein profiles we split the search into chunks
                #I think 32 refers to the number of threads that is running on your computer
  #$file = file that is being accessed in the loop
  #genomad_output_346 = output directory with all 346 MAGs
  #genomad_db = genomad database
  # >> results.out = sending the activity of files through the pipeline to this file to make sure the pipelline is running as it should
  #we could techinically add a --conservative flag here but the default is that sequnces are reuqired to have a plasmid or virus score of at least 0.7 and sequences shorter than 2500bp are required to encorde at least one hallmark gene
    #the conservative flag makes post-classification filters more aggressive preventing sequences without strong support form being classified as plasmid or virus
    #sticking with the default for this case
  genomad end-to-end --cleanup --splits 60 $file genomad_output_final_indiv_assemblies genomad_db >> results.out; 
  done
  
#the files that will be output will be within genomad_output_346. I will need to write a script later on that will combine all corresponding summary files (like phage genes or something) that way we can analyze the data
  
  