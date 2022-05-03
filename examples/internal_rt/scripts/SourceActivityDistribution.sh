Sourcefilelist="$1"
Datafile="$2"
ActivityOutfile="$3"
IsoCenterCoor=""
if [ "$#" -eq 4 ]; then IsoCenterCoor="$4"; fi

declare -a aa

# Origin of the source image
aa=($(cat ${Datafile} | awk '{if($0~/Offset/){print $3" "$4" "$5}}'))
Ox=${aa[0]}; Oy=${aa[1]}; Oz=${aa[2]};
echo "Using Origin = (${Ox}, ${Oy}, ${Oz}) mm"

# Spacing of the source image
aa=($(cat ${Datafile} | awk '{if($0~/ElementSpacing/){print $3" "$4" "$5}}'))
dx=${aa[0]}; dy=${aa[1]}; dz=${aa[2]};
echo "Using Spacing = (${dx}, ${dy}, ${dz}) mm"

# Isocenter source image
Isox="0.0";Isoy="0.0";Isoz="0.0";
if [ ! -z "${IsoCenterCoor}" ]
then
     aa=($(echo ${IsoCenterCoor} | awk '{print $1" "$2" "$3}'))
     Isox=${aa[0]}; Isoy=${aa[1]}; Isoz=${aa[2]};
fi
echo "Using Isocenter = (${Isox}, ${Isoy}, ${Isoz}) mm"

## Origin of the source image
#Ox="281.847"
#Oy="279.547"
#Oz="460.547"

## Spacing of the source image
#dx="4.41806"
#dy="4.41806"
#dz="4.41806"

rm -f ${ActivityOutfile}

rm -f __AUX__

for ifile in ${Sourcefilelist}
do 
        echo -n "Processing ${ifile} ... "

        more ${ifile} | grep "center (image" \
                      | sed 's/voxel center (image source coordinate system):  (/ /g; s/,/ /g; s/)/ /g' \
                      | awk -v Ox="${Ox}" -v Oy="${Oy}" -v Oz="${Oz}" -v Isox="${Isox}" -v Isoy="${Isoy}" -v Isoz="${Isoz}" \
                            -v dx="${dx}" -v dy="${dy}" -v dz="${dz}" '{
                         ix=($1*10-Ox+Isox)/dx
                         iy=($2*10-Oy+Isoy)/dy
                         iz=($3*10-Oz+Isoz)/dz
                         printf("%6.1f %6.1f %6.1f\n",iz,iy,ix) >> "__AUX__"
                        }'
        echo "done"
done
cat __AUX__  | sort -n | uniq -c > __AUX2__
cat __AUX2__ | sort -r -n -k1,1 > ${ActivityOutfile}

rm -f __AUX__ __AUX2__


#tallyimage_spatialSampling:
#                          : voxel index   :   1139233
#                          : source chosen :  ( 33, 68, 69)
#                          : center (image source coordinate system):  (  -13.826,    1.867,  -15.791)
#                          : Position (phantom coordinate system): (   11.460,   12.564,   22.641)
