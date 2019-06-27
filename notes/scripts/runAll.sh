#!/bin/bash

# Arguments
# $1 = dir to work in
# $2 = maf file to use

if [ "$#" -ne 2 ]; then
    echo "Usage: ./runAll.sh [working dir] [MAF alignment]"
    exit 1
fi

# Print a copy of each command to the terminal
# set -x

echo "Writing to $1/log.log"

# Program

# clear out the previous log
rm $1/log.log
touch $1/log.log

# redirect stdout/stderr to log file
exec >> $1/log.log
exec 2>&1

# Split/merge the MAF file into the dir split
echo ""
echo "SPLITTING/MERGING MAF FILE"
echo "" | ./rerunMergeForDirs.sh $1

# cd into the working directory so that we can use the following scripts
cd $1

# Run CMFinder
echo ""
echo "RUNNING CMFINDER"
./../align.sh

# Score with RNAPhylo
echo ""
echo "SCORING CMFINDER RESULTS"
echo "The output for this script is written to the score files output."
./../score.sh ../multiz100tree.newick

# Generate list of RNAPhylo scores for summary
echo ""
echo "GENERATING SCORE SUMMARIES"
echo "The output for this script is written to the listscores output."
./../listScores.sh ./

cd ..

# Generate BED and then bigBed tracks for the found motifs
echo ""
echo "GENERATING BED TRACKS"
mkdir $1/tracks
java -jar trackgenerator.jar -s $1/scores -o $1/tracks/
java -jar reflinegenerator.jar -s $2 -o $1/tracks

# for track hub processing: 
# need to run rerunTracksForDirs separately
# to aggregate the results correctly

