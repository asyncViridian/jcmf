echo "Rerunning trackgenerate etc. on chr1 region"
rm -r q_chr01_5block-10to3000-allspecies/tracks
mkdir q_chr01_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr01_5block-10to3000-allspecies/scores -o q_chr01_5block-10to3000-allspecies/tracks/
for file in q_chr01_5block-10to3000-allspecies/tracks/*.bed
do
  ./generateBigBed.sh $file
done

echo "Rerunning trackgenerate etc. on chr2 region"
rm -r q_chr02_5block-10to3000-allspecies/tracks
mkdir q_chr02_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr02_5block-10to3000-allspecies/scores -o q_chr02_5block-10to3000-allspecies/tracks/
for file in q_chr02_5block-10to3000-allspecies/tracks/*.bed
do
  ./generateBigBed.sh $file
done

echo "Rerunning trackgenerate etc. on chr5 region"
rm -r q_chr05_5block-10to3000-allspecies/tracks
mkdir q_chr05_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr05_5block-10to3000-allspecies/scores -o q_chr05_5block-10to3000-allspecies/tracks/
for file in q_chr05_5block-10to3000-allspecies/tracks/*.bed
do
  ./generateBigBed.sh $file
done

echo "Rerunning trackgenerate etc. on chr12 region"
rm -r q_chr12_5block-10to3000-allspecies/tracks
mkdir q_chr12_5block-10to3000-allspecies/tracks
java -jar trackgenerator.jar -s q_chr12_5block-10to3000-allspecies/scores -o q_chr12_5block-10to3000-allspecies/tracks/
for file in q_chr12_5block-10to3000-allspecies/tracks/*.bed
do
  ./generateBigBed.sh $file
done

