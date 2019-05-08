# Score found alignments with RNAPhylo
echo "Scoring alignments..."
rm -r scores
mkdir scores
cd results
for file in *
do
  RNAPhylo -t ../multiz100tree.newick "$file" > ../scores/"$file".score
done
echo "Scoring done"
