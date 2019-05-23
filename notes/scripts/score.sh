# Score found alignments with RNAPhylo

# Arguments
# $1 = path to tree from working dir

echo "Scoring alignments..."
cp $1 tree.newick
rm -r scores
mkdir scores
cd results
for file in *
do
  RNAPhylo -t ../tree.newick "$file" > ../scores/"$file".score
done
cd ..
rm tree.newick
echo "Scoring done"
