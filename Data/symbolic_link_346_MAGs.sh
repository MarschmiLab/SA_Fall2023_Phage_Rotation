#creating symbolic link of the original binning

for MAGS in 'ls /workdir/lab_data/2014-2015-metaG-muskegon/binning/dereplicated_bins/346_MAGS'
  do
  
  #echo file to make sure that the loop works
  echo $MAGS
  
  #symbolically link the files
  ln -s $MAGS /workdir/sna49/SA_Fall2023_Phage_Rotation/Data

  done
  