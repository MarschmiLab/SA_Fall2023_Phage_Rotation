#set up environment to run genomad
source /programs/miniconda3/bin/activate genomad
genomad download-database .

pwd #make sure that you are in the correct path /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output/

#now time to run the loop to read the files and annonate them!
for file in *.fa.gz; 
  do 
  # Which file am I on? 
  echo $file 
  # input comments here 
  genomad end-to-end --cleanup --splits 32 $file genomad_output_346 genomad_db >>results.out; 
  done
  
#the files that will be output will be within genomad_output_346. I will need to write a script later on that will combine all corresponding summary files (like phage genes or something) that way we can analyze the data
  
  