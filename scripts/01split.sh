# Split the block fasta into individual fasta files
rm -r splitsrc
rm -r split
mkdir splitsrc
mkdir split
echo Splitting FASTA file "$1"
csplit -ksz -f splitsrc/"$1"_ "$1" '/^$/' {*}

# Remove the leading newline from each fasta (csplit aint perfect)
cd splitsrc
for file in *
do
  sed '/./,$!d' "$file" > ../split/"$file"
done
cd ..
rm -r splitsrc
echo Splitting done
