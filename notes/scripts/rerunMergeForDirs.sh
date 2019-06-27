# Re-generate the merged-block MAF files for all working directories passed in
# Example usage: ./rerunMergeForDirs.sh q_chr*

if [ "$#" -eq 0 ]; then
    echo "Usage: ./rerunMergeForDirs [dirs to regen merge for]+"
    exit 1
fi

echo "Please enter addt arguments (by default, runs with no args)"
read args

for dir in "$@"
do
  echo "Rerunning merge in $dir"
  rm -r $dir/split
  mkdir $dir/split
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
  echo ">java -jar blockmerger.jar -s source_mafs/$maf -od $dir/split/ -o m ${args}"
  java -jar blockmerger.jar -s source_mafs/$maf -od $dir/split/ -o m ${args}
  echo
done

