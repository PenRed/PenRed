DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Datafile="$1"
DepthActivitySet1=${2:-"SourceDepthActivity"}
DepthActivitySet2=${3:-"SourceDepthActivity"}
DepthActivityDirSet1=${4:-"PenRed"}
DepthActivityDirSet2=${5:-"Gate"}
TitleSet1=${6:-"PenRed"}
TitleSet2=${7:-"Gate"}

cbcolor=""
if [ "$#" -eq 3 ]; then cbcolor="$3"; fi

declare -a aa

# Dimsize of the source image
aa=($(cat ${Datafile} | awk '{if($0~/DimSize/){print $3" "$4" "$5}}'))
Nx=${aa[0]}; Ny=${aa[1]}; Nz=${aa[2]};
sed 's/__nx__/'${aa[0]}'/;s/__ny__/'${aa[1]}'/;s/__nz__/'${aa[2]}'/' ${DIR}/Comparison-depthdose_template.gnu &> __CACA__

# Number of primary particles generated
#nhisSet1=$(more ${DepthActivityDirSet1}/${DepthActivitySet1}-z.dat |awk 'BEGIN{sum=0}{sum += $2}END{print sum}')
#nhisSet2=$(more ${DepthActivityDirSet2}/${DepthActivitySet2}-z.dat |awk 'BEGIN{sum=0}{sum += $2}END{print sum}')
#sed 's/__nhisSet1__/'${nhisSet1}'/;s/__nhisSet2__/'${nhisSet2}'/' __CACA__ &> __CACA2__

# Depth dose data files
sed 's/__srcdepthSet1__/'${DepthActivitySet1}'/;s/__srcdepthSet2__/'${DepthActivitySet2}'/' __CACA__ &> __CACA2__
sed 's/__srcdepthdirSet1__/'${DepthActivityDirSet1}'/;s/__srcdepthdirSet2__/'${DepthActivityDirSet2}'/' __CACA2__ &> __CACA__
sed 's/__titleSet1__/'${TitleSet1}'/;s/__titleSet2__/'${TitleSet2}'/' __CACA__ &> __CACA2__

gnuplot __CACA2__
rm -f __CACA__ __CACA2__
