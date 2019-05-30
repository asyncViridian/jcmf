# Re-generate the BED/detail tracks for all working directories passed in
# Outputs the resulting bigbed tracks into the temp dir
# Example usage: ./rerunTracksOnAllRegions.sh q_chr*

rm -r temp
mkdir temp

for dir in "$@"
do
  echo "Rerunning trackgenerate etc. in $dir"
  rm -r $dir/tracks
  mkdir $dir/tracks
  java -jar trackgenerator.jar -s $dir/scores -o $dir/tracks
  # Assign the appropriate maf
  if [[ $dir == *"chr01"* ]]; then
    maf="m100_chr1_149844498-149849024.maf"
  elif [[ $dir == *"chr02"* ]]; then
    maf="m100_chr2_218255319-218257366.maf"
  elif [[ $dir == *"chr05"* ]]; then
    maf="m100_chr5_140072857-140108630.maf"
  elif [[ $dir == *"chr12"* ]]; then
    maf="m100_chr12_62602752-62622213.maf"
  fi
  echo "Using MAF file $maf"
  java -jar reflinegenerator.jar -s source_mafs/$maf -o $dir/tracks/
  echo "Appending bed files to sum file"
  cd $dir/tracks/
  for file in *.bed
  do
    cat $file >> ../../temp/$file
  done
  cd ../..
  echo
done

echo "Generating bigBed for sum files"
for file in temp/*.bed
do
  ./generateBigBed.sh $file
done

