#!/bin/bash

# Estimates D_KL(P || Q) (K-L Divergence of P over Q)

# Arguments:
# $1 = dir to work in
# $2 = CM "P"
# $3 = CM "Q"
# $4 = Number of samples to use

# Setup our working environment
mkdir $1/CMs
mkdir $1/emit
mkdir $1/align
cp $2 $1/CMs/
cp $3 $1/CMs/
P=$(basename $2)
Q=$(basename $3)
numSamples=$4
workingDir=$1

# Emit samples
Pemit="${P}.emit.fasta"
./utilities/cmemit -N ${numSamples} -o ${workingDir}/emit/${Pemit} ${workingDir}/CMs/${P}

# Align samples to both CMs
PemitaligntoP="${Pemit}.alignto.${P}.sto"
PemitaligntoQ="${Pemit}.alignto.${Q}.sto"
PemitaligntoPdata="${Pemit}.alignto.${P}.sto.data"
PemitaligntoQdata="${Pemit}.alignto.${Q}.sto.data"
./utilities/cmalign -o ${workingDir}/align/${PemitaligntoP} ${workingDir}/CMs/${P} ${workingDir}/emit/${Pemit} > ${workingDir}/align/${PemitaligntoPdata}
./utilities/cmalign -o ${workingDir}/align/${PemitaligntoQ} ${workingDir}/CMs/${Q} ${workingDir}/emit/${Pemit} > ${workingDir}/align/${PemitaligntoQdata}

# Read the bitscores of all alignments and calculate score
# read a file line by line:
while IFS= read -r -u3 line; do
    if [[ ${line} != \#* ]]; then
        # Ignore the comment lines in data output
        echo "File: ${line}"
    fi
done 3< "${workingDir}/align/${PemitaligntoPdata}"

echo "TODO"

