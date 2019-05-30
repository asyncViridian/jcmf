# Score found alignments with RNAPhylo

# Arguments
# $1 = path to tree from working dir

echo "Scoring alignments..."
rm -r scores
mkdir scores
cd results
for file in *
do
  RNAPhylo -t ../$1 --partition "$file" > ../scores/"$file".score
done
cd ..
echo "Scoring done"
