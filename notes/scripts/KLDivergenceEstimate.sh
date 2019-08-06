#!/bin/bash

# Estimates D_KL(P || Q) (K-L Divergence of P over Q)

# Arguments:
# $1 = dir to work in
# $2 = CM "P"
# $3 = CM "Q"
# $4 = Number of samples to use

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

# Variables for the simple average method:
# using only emits from P
pAlignSum=0
# using emits from both P and Q
pqAlignSum=0

# Variables for the 2^Pi (2^a) adjusted method:
# using only emits from P
aAdjustPAlignNum=0
aAdjustPAlignDen=0
# using emits from both P and Q
aAdjustPQAlignNum=0
aAdjustPQAlignDen=0

# read a file line by line:
while IFS= read -r -u3 line; do
    if [[ ${line} != \#* ]]; then
        # Ignore the comment lines in data output
        Plinearr=(${line})
        Qlinearr=($(grep ${Plinearr[1]} ${workingDir}/align/${PemitaligntoQdata}))
        # Bit score is arr[6] in the whitespace-sv
        echo "P ${Plinearr[6]} Q ${Qlinearr[6]}"
        # To take the power of x^n : e($n*l($x))
        PQ=$(echo "(${Plinearr[6]})-(${Qlinearr[6]})" | bc -l)
        #echo "Elem score: ${PQ}"
        pAlignSum=$(echo "(${pAlignSum})+(${PQ})" | bc -l)
        pqAlignSum=$(echo "(${pqAlignSum})+(${PQ})" | bc -l)
        aAdjustPAlignNum=$(echo "${aAdjustPAlignNum}+((e((${Plinearr[6]})*l(2)))*(${PQ}))" | bc -l)
        aAdjustPAlignDen=$(echo "${aAdjustPAlignDen}+((e((${Plinearr[6]})*l(2))))" | bc -l)
        aAdjustPQAlignNum=$(echo "${aAdjustPQAlignNum}+((e((${Plinearr[6]})*l(2)))*(${PQ}))" | bc -l)
        aAdjustPQAlignDen=$(echo "${aAdjustPQAlignDen}+((e((${Plinearr[6]})*l(2))))" | bc -l)
    fi
done 3< "${workingDir}/align/${PemitaligntoPdata}"
while IFS= read -r -u3 line; do
    if [[ ${line} != \#* ]]; then
        # Ignore the comment lines in data output
        Plinearr=(${line})
        Qlinearr=($(grep ${Plinearr[1]} ${workingDir}/align/${QemitaligntoQdata}))
        # Bit score is arr[6] in the whitespace-sv
        echo "P ${Plinearr[6]} Q ${Qlinearr[6]}"
        # To take the power of x^n : e($n*l($x))
        PQ=$(echo "(${Plinearr[6]})-(${Qlinearr[6]})" | bc -l)
        #echo "Elem score: ${PQ}"
        pqAlignSum=$(echo "(${pqAlignSum})+(${PQ})" | bc -l)
        aAdjustPQAlignNum=$(echo "${aAdjustPQAlignNum}+((e((${Plinearr[6]})*l(2)))*(${PQ}))" | bc -l)
        aAdjustPQAlignDen=$(echo "${aAdjustPQAlignDen}+((e((${Plinearr[6]})*l(2))))" | bc -l)
    fi
done 3< "${workingDir}/align/${QemitaligntoPdata}"


# Average scores (using sum and sample count)
#echo "sum ${pAlignSum} num ${numSamples}"
echo "uncombine reg $(echo "${pAlignSum}/${numSamples}" | bc -l)"
echo "uncombine cor $(echo "${aAdjustPAlignNum}/${aAdjustPAlignDen}" | bc -l)"
echo "  combine reg $(echo "${pqAlignSum}/(${numSamples}*2)" | bc -l)"
echo "  combine cor $(echo "${aAdjustPQAlignNum}/(${aAdjustPQAlignDen})" | bc -l)"

