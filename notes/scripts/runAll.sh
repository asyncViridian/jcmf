#!/bin/bash

# Arguments
# $1 = dir to work in
# $2 = maf file to use

exec > $1/log 2>&1
set -x

if [ "$#" -ne 2 ]; then
    echo "Usage: ./runAll.sh [working dir] [MAF alignment]"
    exit 1
fi

# Program

# clear out the previous log
rm $1/log

# make the split directory to put fasta chunks in
rm -r $1/split
mkdir $1/split

# Split/merge the MAF file into the dir split
echo ""
echo "SPLITTING/MERGING MAF FILE"
java -jar blockmerger.jar -s $2 -od $1/split -o m

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

# Generate BED and then bigBed tracks for the found motifs
echo ""
echo "GENERATING BED TRACKS"
cd ..
mkdir $1/tracks
java -jar trackgenerator.jar -s $1/scores -o $1/tracks/
for file in $1/tracks/*
do
  ./generateBigBed.sh file
done

# TODO add track hub processing?

