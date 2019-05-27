rm -r temp
mkdir temp

echo "Rerunning trackgenerate etc. on chr1 region"
rm -r q_chr01_5block-10to3000-allspecies/tracks
mkdir q_chr01_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr01_5block-10to3000-allspecies/scores -o q_chr01_5block-10to3000-allspecies/tracks/
java -jar reflinegenerator.jar -s source_mafs/m100_chr1_149844498-149849024.maf -o q_chr01_5block-10to3000-allspecies/tracks/
echo "Appending bed files to sum file"
cd q_chr01_5block-10to3000-allspecies/tracks/
for file in *.bed
do
  cat $file >> ../../temp/$file
done
cd ../..

echo "Rerunning trackgenerate etc. on chr2 region"
rm -r q_chr02_5block-10to3000-allspecies/tracks
mkdir q_chr02_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr02_5block-10to3000-allspecies/scores -o q_chr02_5block-10to3000-allspecies/tracks/
java -jar reflinegenerator.jar -s source_mafs/m100_chr2_218255319-218257366.maf -o q_chr02_5block-10to3000-allspecies/tracks/
echo "Appending bed files to sum file"
cd q_chr02_5block-10to3000-allspecies/tracks/
for file in *.bed
do
  cat $file >> ../../temp/$file
done
cd ../..

echo "Rerunning trackgenerate etc. on chr5 region"
rm -r q_chr05_5block-10to3000-allspecies/tracks
mkdir q_chr05_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr05_5block-10to3000-allspecies/scores -o q_chr05_5block-10to3000-allspecies/tracks/
java -jar reflinegenerator.jar -s source_mafs/m100_chr5_140072857-140108630.maf -o q_chr05_5block-10to3000-allspecies/tracks/
echo "Appending bed files to sum file"
cd q_chr05_5block-10to3000-allspecies/tracks/
for file in *.bed
do
  cat $file >> ../../temp/$file
done

echo "Rerunning trackgenerate etc. on chr12 region"
rm -r q_chr12_5block-10to3000-allspecies/tracks
mkdir q_chr12_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr12_5block-10to3000-allspecies/scores -o q_chr12_5block-10to3000-allspecies/tracks/
java -jar reflinegenerator.jar -s source_mafs/m100_chr12_62602752-62622213.maf -o q_chr12_5block-10to3000-allspecies/tracks/
echo "Appending bed files to sum file"
cd q_chr02_5block-10to3000-allspecies/tracks/
for file in *.bed
do
  cat $file >> ../../temp/$file
done
cd ../..

echo "Generating bigBed for sum files"
for file in temp/*.bed
do
  ./generateBigBed.sh $file
done

