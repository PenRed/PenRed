DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Datafile="$1"
zslices="$2"
title1="$3"
title2="$4"
outname1="$5"
outname2="$6"
cbcolor=""
if [ "$#" -eq 7 ]; then cbcolor="$7"; fi

declare -a aa

# Dimsize of the source image
aa=($(cat ${Datafile} | awk '{if($0~/DimSize/){print $3" "$4" "$5}}'))
Nx=${aa[0]}; Ny=${aa[1]}; Nz=${aa[2]};
sed 's/__nx__/'${aa[0]}'/;s/__ny__/'${aa[1]}'/' ${DIR}/plotImage_template.gplot &> __CACA-11__
sed 's/__nx__/'${aa[0]}'/;s/__ny__/'${aa[1]}'/' ${DIR}/plotImage-png_template.gplot &> __CACA-21__

# Z slices to show
aa=($(echo ${zslices} | awk '{print $1" "$2" "$3" "$4" "$5" "$6}'))
sed 's/__pl1__/'${aa[0]}'/;s/__pl2__/'${aa[1]}'/;s/__pl3__/'${aa[2]}'/;s/__pl4__/'${aa[3]}'/;s/__pl5__/'${aa[4]}'/;s/__pl6__/'${aa[5]}'/' __CACA-11__ &> __CACA-12__    
sed 's/__pl1__/'${aa[0]}'/;s/__pl2__/'${aa[1]}'/;s/__pl3__/'${aa[2]}'/;s/__pl4__/'${aa[3]}'/;s/__pl5__/'${aa[4]}'/;s/__pl6__/'${aa[5]}'/' __CACA-21__ &> __CACA-22__    
echo "Using zslices = (${aa[@]})"

# Cbcolors to show
cb1="*";cb2="*";cb3="*";cb4="*";cb5="*";cb6="*";
if [ ! -z "${cbcolor}" ]
then
     aa=($(echo ${cbcolor} | awk '{print $1" "$2" "$3" "$4" "$5" "$6}'))
     cb1=${aa[0]}; cb2=${aa[1]}; cb3=${aa[2]}; cb4=${aa[3]}; cb5=${aa[4]}; cb6=${aa[5]};
fi
sed 's/__cb1__/'${cb1}'/;s/__cb2__/'${cb2}'/;s/__cb3__/'${cb3}'/;s/__cb4__/'${cb4}'/;s/__cb5__/'${cb5}'/;s/__cb6__/'${cb6}'/' __CACA-12__ &> __CACA-11__    
sed 's/__cb1__/'${cb1}'/;s/__cb2__/'${cb2}'/;s/__cb3__/'${cb3}'/;s/__cb4__/'${cb4}'/;s/__cb5__/'${cb5}'/;s/__cb6__/'${cb6}'/' __CACA-22__ &> __CACA-21__    
echo "Using cbcolors = (${cb1}, ${cb2}, ${cb3}, ${cb4}, ${cb5}, ${cb6})"

# Titles
sed 's/__title1__/'"${title1}"'/g;s/__title2__/'"${title2}"'/g' __CACA-11__ &> __CACA-12__
sed 's/__title1__/'"${title1}"'/g;s/__title2__/'"${title2}"'/g' __CACA-21__ &> __CACA-22__
echo "Using title 1 = ${title1}"
echo "Using title 2 = ${title2}"

# Output files
sed 's/__outnameSet1__/'"${outname1}"'/g;s/__outnameSet2__/'"${outname2}"'/g' __CACA-12__ &> __CACA-11__
sed 's/__outnameSet1__/'"${outname1}"'/g;s/__outnameSet2__/'"${outname2}"'/g' __CACA-22__ &> __CACA-21__

gnuplot __CACA-11__
gnuplot __CACA-21__
rm -f __CACA-11__ __CACA-21__ __CACA-12__ __CACA-22__
