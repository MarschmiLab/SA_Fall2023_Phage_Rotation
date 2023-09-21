#set up environment to run genomad
source /programs/miniconda3/bin/activate genomad
genomad download-database .

#make sure that you are in the correct path /local/workdir/sna49/SA_Fall2023_Phage_Rotation/Data/346_MAGs/Output
for file in *.fa.gz; 
  do 
  # Which file am I on? 
  echo $file 
  # input comments here 
  genomad end-to-end --cleanup --splits 32 $file genomad_output genomad_db >>results.out; 
  done
  