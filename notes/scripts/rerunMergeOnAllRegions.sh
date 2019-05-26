echo "Rerunning merge on chr1 region"
rm -r q_chr01_5block-10to3000-allspecies/split
mkdir q_chr01_5block-10to3000-allspecies/split
java -jar blockmerger.jar -s source_mafs/m100_chr1_149844498-149849024.maf -od q_chr01_5block-10to3000-allspecies/split/ -o m

echo "Rerunning merge on chr2 region"
rm -r q_chr02_5block-10to3000-allspecies/split
mkdir q_chr02_5block-10to3000-allspecies/split
java -jar blockmerger.jar -s source_mafs/m100_chr2_218255319-218257366.maf -od q_chr02_5block-10to3000-allspecies/split/ -o m

echo "Rerunning merge on chr5 region"
rm -r q_chr05_5block-10to3000-allspecies/split
mkdir q_chr05_5block-10to3000-allspecies/split
java -jar blockmerger.jar -s source_mafs/m100_chr5_140072857-140108630.maf -od q_chr05_5block-10to3000-allspecies/split/ -o m

echo "Rerunning merge on chr12 region"
rm -r q_chr12_5block-10to3000-allspecies/split
mkdir q_chr12_5block-10to3000-allspecies/split
java -jar blockmerger.jar -s source_mafs/m100_chr12_62602752-62622213.maf -od q_chr12_5block-10to3000-allspecies/split/ -o m


