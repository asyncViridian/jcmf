# Generate bigBed files from the given bed files

# Arguments
# $1 = bed file to convert

echo "Sorting $1 ..."
./utilities/bedSort $1 $1
echo "Converting $1 ..."
./utilities/bedToBigBed $1 http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes ${1%.bed}.bb
echo "Converted $1 -> ${1%.bed}.bb"

