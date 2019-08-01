#!/bin/bash

# Arguments
# $1 = dir to work in
# $2 = maf file to use
# $3 = settings file to use

if [ "$#" -ne 3 ]; then
    echo "Usage: ./runAll.sh [working dir] [MAF alignment] [settings file]"
    exit 1
fi

# Print a copy of each command to the terminal
# set -x

# Program

# clear out the previous log
rm $1/log.log
touch $1/log.log

# redirect stdout/stderr to log file
exec >> $1/log.log
exec 2>&1

# Get variables from settings file
echo "-=- GETTING SETTINGS VARIABLES"
source $3
echo "mergeArgs=\"${mergeArgs}\""
echo "maxCandPerSeq=\"${maxCandPerSeq}\""
echo "fracSeqWMotif=\"${fracSeqWMotif}\""
echo "maxNumSingleS=\"${maxNumSingleS}\""
echo "maxNumDoubleS=\"${maxNumDoubleS}\""
echo "minSpanSingle=\"${minSpanSingle}\""
echo "maxSpanSingle=\"${maxSpanSingle}\""
echo "minSpanDouble=\"${minSpanDouble}\""
echo "maxSpanDouble=\"${maxSpanDouble}\""
echo "minCandLength=\"${minCandLength}\""
echo "maxCandLength=\"${maxCandLength}\""

# Split/merge the MAF file into the dir split
echo ""
echo "-=- SPLITTING/MERGING MAF FILE"
rm -r $1/split*
mkdir $1/split
java -jar blockmerger.jar -s $2 -od $1/split -o m -ot maf ${mergeArgs}

# DIFFERENT FROM ORIGINAL:
# SHUFFLE THE MAF MERGEDBLOCKS
mkdir $1/splitshuffle
cd $1/splitshuffle
for file in ../split/*.maf
do
    ./../../utilities/multiperm --num=10 ${file}
done
for file in ../split/*.png
do
    mv $file $(basename ${file})
done
cd ../..
# MOVE THE SHUFFLED BLOCKS INTO SPLIT
rm -r $1/split
mv $1/splitshuffle $1/split
# CONVERT THE SHUFFLED MAFS INTO FASTA
mkdir $1/splitfasta
for file in $1/split/*.maf
do
    filename=$(basename ${file})
    cat ${file} | awk '/^s/{print ">" $2 "\n" $7}' | sed 's/-//g' >> $1/splitfasta/${filename%.maf}.fasta
done
for file in $1/split/*.png
do
    mv $file $1/splitfasta/$(basename ${file})
done
rm -r $1/split
mv $1/splitfasta $1/split

# cd into the working directory so that we can use the following scripts
cd $1

# Run CMFinder
echo ""
echo "-=- RUNNING CMFINDER"
./../align.sh ../$3

# Score with RNAPhylo
echo ""
echo "-=- SCORING CMFINDER RESULTS"
echo "The output for this script is written to the score files output."
./../score.sh ../multiz100tree.newick

# Generate list of RNAPhylo scores for summary
echo ""
echo "-=- GENERATING SCORE SUMMARIES"
echo "The output for this script is written to the listscores output."
./../listScores.sh ./

cd ..

# Generate BED and then bigBed tracks for the found motifs
echo ""
echo "-=- GENERATING BED TRACKS"
mkdir $1/tracks
java -jar trackgenerator.jar -s $1/scores -o $1/tracks/
java -jar reflinegenerator.jar -s $2 -o $1/tracks

# for track hub processing: 
# need to run rerunTracksForDirs separately
# to aggregate the results correctly

