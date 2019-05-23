#!/bin/sh

# Arguments
# $1 = dir to work in
# $2 = maf file to use

if [ "$#" -ne 2 ]; then
    echo "Usage: ./runAll.sh [working dir] [MAF alignment]"
    exit 1
fi

# Program

# make the split directory to put fasta chunks in
rm -r $1/split
mkdir $1/split

# Split/merge the MAF file into the dir split
java -jar blockmerger.jar -s $2 -od $1/split -o m >> $1/log

# cd into the working directory so that we can use the following scripts
cd $1

# Run CMFinder
./../align.sh

# Score with RNAPhylo
./../score.sh

# Generate list of RNAPhylo scores for summary
./../listScores.sh ./

# Generate BED tracks for the found motifs
cd ..
mkdir $1/tracks
java -jar trackgenerator.jar $1/scores -o $1/tracks/multi.bed # TODO revise?

# TODO add bigBed processing
