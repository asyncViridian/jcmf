#!/bin/bash

# Estimates D_KL(P || Q) (K-L Divergence of P over Q)
# OUTPUTS 6 SEPARATE SCORES, EACH COMPUTED USING A SEPARATE METHOD

# Arguments:
# $1 = dir to work in
# $2 = CM "P"
# $3 = CM "Q"
# $4 = Number of samples to use emitted from each of P, Q
#      (So if using a combined method, the actual number of samples is 2x)

# Setup our working environment
mkdir -p $1/CMs
mkdir -p $1/emit
mkdir -p $1/align
cp $2 $1/CMs/
cp $3 $1/CMs/
P=$(basename $2)
Q=$(basename $3)
numSamples=$4
workingDir=$1

# Get stats about the CMs
tmp=($(cat ${workingDir}/CMs/${P} | grep "CLEN"))
Plen=${tmp[1]}
tmp=($(cat ${workingDir}/CMs/${Q} | grep "CLEN"))
Qlen=${tmp[1]}

# Emit samples
Pemit="${P}.emit.fasta"
Qemit="${Q}.emit.fasta"
./utilities/cmemit -N ${numSamples} -o ${workingDir}/emit/${Pemit} ${workingDir}/CMs/${P}
./utilities/cmemit -N ${numSamples} -o ${workingDir}/emit/${Qemit} ${workingDir}/CMs/${Q}

# Align samples to both CMs
PemitaligntoP="${Pemit}.alignto.${P}.sto"
PemitaligntoQ="${Pemit}.alignto.${Q}.sto"
PemitaligntoPdata="${Pemit}.alignto.${P}.sto.data"
PemitaligntoQdata="${Pemit}.alignto.${Q}.sto.data"
./utilities/cmalign -o ${workingDir}/align/${PemitaligntoP} ${workingDir}/CMs/${P} ${workingDir}/emit/${Pemit} > ${workingDir}/align/${PemitaligntoPdata}
./utilities/cmalign -o ${workingDir}/align/${PemitaligntoQ} ${workingDir}/CMs/${Q} ${workingDir}/emit/${Pemit} > ${workingDir}/align/${PemitaligntoQdata}
QemitaligntoP="${Qemit}.alignto.${P}.sto"
QemitaligntoQ="${Qemit}.alignto.${Q}.sto"
QemitaligntoPdata="${Qemit}.alignto.${P}.sto.data"
QemitaligntoQdata="${Qemit}.alignto.${Q}.sto.data"
./utilities/cmalign -o ${workingDir}/align/${QemitaligntoP} ${workingDir}/CMs/${P} ${workingDir}/emit/${Qemit} > ${workingDir}/align/${QemitaligntoPdata}
./utilities/cmalign -o ${workingDir}/align/${QemitaligntoQ} ${workingDir}/CMs/${Q} ${workingDir}/emit/${Qemit} > ${workingDir}/align/${QemitaligntoQdata}

# Read the bitscores of all alignments and calculate score

# Variables for the simple average method (exaggerates differences)
# using only emits from P
pAlignSum=0
# using emits from both P and Q
pqAlignSum=0

# Variables for the 2^Pi (2^a) adjusted method (should correct somewhat for emit distribution)
# using only emits from P
aAdjustPAlignNum=0
aAdjustPAlignDen=0
# using emits from both P and Q
aAdjustPQAlignNum=0
aAdjustPQAlignDen=0

# Variables for the 2^Pi/N (2^a/N) adjusted method (should be more stable over nulls)
# using only emits from P
nAdjustPAlignNum=0
nAdjustPAlignDen=0
# using emits from both P and Q
nAdjustPQAlignNum=0
nAdjustPQAlignDen=0

# read a file line by line:
while IFS= read -r -u3 line; do
    if [[ ${line} != \#* ]]; then
        # Ignore the comment lines in data output
        Plinearr=(${line})
        Qlinearr=($(grep ${Plinearr[1]} ${workingDir}/align/${PemitaligntoQdata}))
        # Emit seq length is arr[2] in the whitespace-sv
        emitLen=${Plinearr[2]}
        # Bit score is arr[6] in the whitespace-sv
        Pbitsc=${Plinearr[6]}
        Qbitsc=${Qlinearr[6]}
        #echo "P ${Pbitsc} Q ${Qbitsc}"
        PQ=$(echo "(${Pbitsc})-(${Qbitsc})" | bc -l)
        #echo "Elem score: ${PQ}"
        pAlignSum=$(echo "(${pAlignSum})+(${PQ})" | bc -l)
        pqAlignSum=$(echo "(${pqAlignSum})+(${PQ})" | bc -l)
        aterm=$(echo "e((${Pbitsc})*l(2))" | bc -l)
        aAdjustPAlignNum=$(echo "${aAdjustPAlignNum}+((${aterm})*(${PQ}))" | bc -l)
        aAdjustPAlignDen=$(echo "${aAdjustPAlignDen}+((${aterm}))" | bc -l)
        aAdjustPQAlignNum=$(echo "${aAdjustPQAlignNum}+((${aterm})*(${PQ}))" | bc -l)
        aAdjustPQAlignDen=$(echo "${aAdjustPQAlignDen}+((${aterm}))" | bc -l)
        nterm=$(echo "((1/4)^(${emitLen}-${Plen}))*(e(-1*(${emitLen}/${Plen})))" | bc -l)
        nAdjustPAlignNum=$(echo "${nAdjustPAlignNum}+((${aterm})*(${nterm})*(${PQ}))" | bc -l)
        nAdjustPAlignDen=$(echo "${nAdjustPAlignDen}+((${aterm})*(${nterm}))" | bc -l)
        nAdjustPQAlignNum=$(echo "${nAdjustPQAlignNum}+((${aterm})*(${nterm})*(${PQ}))" | bc -l)
        nAdjustPQAlignDen=$(echo "${nAdjustPQAlignDen}+((${aterm})*(${nterm}))" | bc -l)
        # To take the power of x^n : e($n*l($x))
    fi
done 3< "${workingDir}/align/${PemitaligntoPdata}"
while IFS= read -r -u3 line; do
    if [[ ${line} != \#* ]]; then
        # Ignore the comment lines in data output
        Plinearr=(${line})
        Qlinearr=($(grep ${Plinearr[1]} ${workingDir}/align/${QemitaligntoQdata}))
        # Emit seq length is arr[2] in the whitespace-sv
        emitLen=${Qlinearr[2]}
        # Bit score is arr[6] in the whitespace-sv
        Pbitsc=${Plinearr[6]}
        Qbitsc=${Qlinearr[6]}
        #echo "P ${Pbitsc} Q ${Qbitsc}"
        PQ=$(echo "(${Pbitsc})-(${Qbitsc})" | bc -l)
        #echo "Elem score: ${PQ}"
        pqAlignSum=$(echo "(${pqAlignSum})+(${PQ})" | bc -l)
        aterm=$(echo "e((${Pbitsc})*l(2))" | bc -l)
        aAdjustPQAlignNum=$(echo "${aAdjustPQAlignNum}+((${aterm})*(${PQ}))" | bc -l)
        aAdjustPQAlignDen=$(echo "${aAdjustPQAlignDen}+((${aterm}))" | bc -l)
        nterm=$(echo "((1/4)^(${emitLen}-${Qlen}))*(e(-1*(${emitLen}/${Qlen})))" | bc -l)
        nAdjustPQAlignNum=$(echo "${nAdjustPQAlignNum}+((${aterm})*(${nterm})*(${PQ}))" | bc -l)
        nAdjustPQAlignDen=$(echo "${nAdjustPQAlignDen}+((${aterm})*(${nterm}))" | bc -l)
    fi
done 3< "${workingDir}/align/${QemitaligntoPdata}"


# Average scores (using sum and sample count)
# uncombined (P-only) direct avg
echo "$(echo "${pAlignSum}/${numSamples}" | bc -l)"
# combined (P+Q) direct avg
echo "$(echo "${pqAlignSum}/(${numSamples}*2)" | bc -l)"
# uncombined (P-only) corrected with 2^Pi
echo "$(echo "${aAdjustPAlignNum}/${aAdjustPAlignDen}" | bc -l)"
# combined (P+Q) corrected with 2^Pi
echo "$(echo "${aAdjustPQAlignNum}/(${aAdjustPQAlignDen})" | bc -l)"
# uncombined (P-only) corrected with 2^Pi/N
echo "$(echo "${nAdjustPAlignNum}/${nAdjustPAlignDen}" | bc -l)"
# combined (P+Q) corrected with 2^Pi/N
echo "$(echo "${nAdjustPQAlignNum}/(${nAdjustPQAlignDen})" | bc -l)"
