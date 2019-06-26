# arg1 = directory to get all images from
# arg2 = directory to put PNGs in

# get an id (the same as the pid of this script) to match the copyphoto script this generates
id=$$

find $1 -regex '.*\(jpg\|jpeg\|png\|gif\)' \! -path './$2/*' -exec echo cp -t $2 {} + > copyphoto${id}.sh
chmod +x copyphoto${id}.sh
./copyphoto${id}.sh
rm copyphoto${id}.sh
