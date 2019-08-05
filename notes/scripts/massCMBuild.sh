#!/bin/bash

# Converts all of the CMfinder result files (scored MSAs) 
# in a given directory into CMs

# Arguments:
# $1 = dir where all the CMfinder outputs are (no particular extension)
# $2 = dir to write all the converted CMs to (with .cm extension)

mkdir -p $2/info
for file in $1/*
do
    cmfile=$(basename ${file}).cm
    ./utilities/cmbuild -F $2/${cmfile} ${file} > $2/info/${cmfile}.info
done

