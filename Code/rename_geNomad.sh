#compressing all fasta files to run with geNomad

#all the files from the metagenomic data end in .fa and we are writing a script that will add .gz at the end of all the files

#make sure that you are in the correct directory, we should be in the /346_MAGs/Output directory
for file in *.fa
do
  mv "$file" "${file%.fa}.fa.gz"
#gzip is needed to compress all the files
  gzip "${file%.fa}.fa.gz"
  echo "renamed files to .gz"
done

echo "all done!"
