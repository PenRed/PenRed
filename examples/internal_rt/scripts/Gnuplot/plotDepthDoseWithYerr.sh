DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Datafile="$1"
DepthDosePenR="$2"
DepthDoseGate="$3"

cbcolor=""
if [ "$#" -eq 3 ]; then cbcolor="$3"; fi

declare -a aa

# Dimsize of the source image
aa=($(cat ${Datafile} | awk '{if($0~/DimSize/){print $3" "$4" "$5}}'))
Nx=${aa[0]}; Ny=${aa[1]}; Nz=${aa[2]};
sed 's/__nx__/'${aa[0]}'/;s/__ny__/'${aa[1]}'/;s/__nz__/'${aa[2]}'/' ${DIR}/Comparison-depthdose_wyerr_template.gnu &> __CACA__

# Number of primary particles generated
#nhisPenR=$(more PenRed/${DepthDosePenR}-z.dat |awk 'BEGIN{sum=0}{sum += $2}END{print sum}')
#nhisGate=$(more Gate/${DepthDoseGate}-z.dat |awk 'BEGIN{sum=0}{sum += $2}END{print sum}')
#sed 's/__nhisPenR__/'${nhisPenR}'/;s/__nhisGate__/'${nhisGate}'/' __CACA__ &> __CACA2__

# Depth dose data files
sed 's/__srcdepthPenR__/'${DepthDosePenR}'/;s/__srcdepthGate__/'${DepthDoseGate}'/' __CACA__ &> __CACA2__

gnuplot __CACA2__
rm -f __CACA__ __CACA2__
