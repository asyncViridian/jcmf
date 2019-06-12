# arg1 = directory to get all images from
# arg2 = directory to put PNGs in
find $1 -regex '.*\(jpg\|jpeg\|png\|gif\)' \! -path './$2/*' -exec echo cp -t $2 {} + > copyphoto.sh
chmod +x copyphoto.sh
./copyphoto.sh
rm copyphoto.sh
