#!/bin/bash

# Arguments:
# $1 = dir where all the CMs are (with file extension .cm)
# $2 = Number of samples to use for divergence calculation

# Generates a whitespace-delimited matrix of K-L divergences
# D_KL(P || Q)
# for all CMs within a given directory. Outputs to stdout.

# Output format: Each entry in the matrix = D_KL(P || Q)
# #    cm1  cm2  ...  cmN
# cm1  1||1 1||2 ...  1||N
# cm2  2||1 2||2 ...  2||N
# ...  ...  ...  ...  ... 
# cmN  N||1 N||2 ...  N||N

CMs=$1/*.cm
tempDir=crossscore$$
numSamples=$2

# Print the header line
printf "#header \t"
for fileQ in ${CMs}
do
    printf "$(basename ${fileQ}) \t"
done
printf "\n"

# Create temp directory to work in
mkdir -p ${tempDir}
# Print the actual table contents
for fileP in ${CMs}
do
    printf "$(basename ${fileP}) \t"
    for fileQ in ${CMs}
    do
        # Generate K-L score
        PQdiv=$(./KLDivergenceEstimate.sh ${tempDir}/ ${fileP} ${fileQ} ${numSamples})
        printf "${PQdiv} \t"
    done
    printf "\n"
done

rm -r ${tempDir}
