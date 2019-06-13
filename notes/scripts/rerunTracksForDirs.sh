# Re-generate the BED/detail tracks for all working directories passed in
# Outputs the resulting bigbed tracks into the temp dir
# Example usage: ./rerunTracksOnAllRegions.sh q_chr*

echo "Please enter unique identifier to prefix outputs with"
read prefix

echo "Re-generating temp directory"
rm -r temp_$prefix
mkdir temp_$prefix

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
    cat $file >> ../../temp_$prefix/${prefix}_$file
  done
  cd ../..
  echo
done

echo "Generating bigBed for sum files"
for file in temp_$prefix/*.bed
do
  ./generateBigBed.sh $file
done

echo "Creating source references files"
# Copying the source files to a temp dir
mkdir temp_$prefix/tempsrc
for dir in "$@"
do
  echo "Copying over all source files for $dir"
  # Assign the appropriate chr postfix
  if [[ $dir == *"chr01"* ]]; then
    chr="chr1"
  elif [[ $dir == *"chr02"* ]]; then
    chr="chr2"
  elif [[ $dir == *"chr05"* ]]; then
    chr="chr5"
  elif [[ $dir == *"chr12"* ]]; then
    chr="chr12"
  fi
  cd $dir/scores
  for file in *.score
  do
    cp $file ../../temp_$prefix/tempsrc/${file}_${chr}.txt
  done
  cd ../..
done
# Grabbing only the source files that we actually want
mkdir temp_$prefix/${prefix}_src
cd temp_$prefix
for file in *.bed
do
  # Don't copy over source files for the reference sequence ._.
  if [[ ! $file =~ "ref" ]]; then
    # Split each line of the bed files into array
    while IFS=$'\t' read -r -a lineArray
    do
      filename="${lineArray[3]}_${lineArray[0]}.txt"
      echo "Saving to src: $filename"
      cp tempsrc/$filename ${prefix}_src/$filename
    done < $file
  fi
done
cd ..
rm -r temp_$prefix/tempsrc

# Getting the PNG files for each source dir
mkdir temp_$prefix/graphics
for dir in "$@"
do
  echo "Copying PNGs out of $dir"
  mkdir temp_$prefix/graphics/$dir
  ./collectPNGs.sh $dir temp_$prefix/graphics/$dir
done
