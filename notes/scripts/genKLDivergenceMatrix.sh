#!/bin/bash

# Generates a whitespace-delimited matrix of symmetrical K-L divergences
# D_KL(P || Q) + D_KL(Q || P)
# for all CMs within a given directory. Outputs to stdout.

# Arguments:
# $1 = dir where all the CMs are (with file extension .cm)
# $2 = Number of samples to use for divergence calculation

for fileP in $1/*.cm
do
    for fileQ in $1/*.cm
    do
        echo "CMfiles ${fileP} ${fileQ}"
    done
done

