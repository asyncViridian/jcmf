# Run CMFinder (log output in log file, wipe log file on each run)

# Flag vars
maxCandPerSeq=60 # pl-def 40
fracSeqWMotif=0.70 # pl-def 0.8
maxNumSingleS=20 # pl-def 5
maxNumDoubleS=20 # pl-def 5
minSpanSingle=20 # pl-def 30
maxSpanSingle=150 # pl-def 100
minSpanDouble=20 # pl-def 40
maxSpanDouble=150 # pl-def 100
minCandLength=30 # Unused as an input flag for cmfinder.pl
maxCandLength=150 # Unused as an input flag for cmfinder.pl

touch log
rm log
touch log
rm -r times
mkdir times
cd times
mkdir split
cd ..

echo alignment running
for file in split/*.fasta; do
  # Prevent the loop from trying to run on "*.fasta" if no files exist
  [ -f "$file" ] || break
  # Run CMFinder
  cmfinder04.pl \
                -c $maxCandPerSeq \
                -f $fracSeqWMotif \
                -s1 $maxNumSingleS \
                -s2 $maxNumDoubleS \
                -minspan1 $minSpanSingle \
                -maxspan1 $maxSpanSingle \
                -minspan2 $minSpanDouble \
                -maxspan2 $maxSpanDouble \
                -combine \
                -amaa \
                -allCpus \
                -candsParallel \
                -saveTimer "times/$file" \
                "$file" >> log 2>&1
                #-likeold \ # Why does this return 0 results?
                #-m $minCandLength \ # This is not a real arg...
                #-M $maxCandLength \ # This is not a real arg...
done
echo alignment done

# Grab successfully found motifs from CMfinder and move them to results/*
rm -r results
mkdir results
mv split/*motif* results/

# Remove other side products of alignment so we don't deal with extra files
rm -r other
mkdir other
mv split/*file-list* other/
mv split/*cand* other/
mv split/*align* other/

# Remove temp and dnd files from the results too
cd results
mv *.temp* ../other/
mv *.dnd* ../other/
cd ..

# Clean up time output
cd times
mv split/* ./
rmdir split
cd ..
