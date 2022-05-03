Activityfile="$1"
Datafile="$2"
Headeroutfile="$3"

# Dimsize of the source image
aa=($(cat ${Datafile} | awk '{if($0~/DimSize/){print $3" "$4" "$5}}'))
Nx=${aa[0]}; Ny=${aa[1]}; Nz=${aa[2]};

#Nx=128
#Ny=128
#Nz=209

for i in `seq 0 $(((${Nx}-1)))`
do 
   more ${Activityfile} | awk -v ind="$i" 'BEGIN{sum=0}{if($4==ind){sum += $1}}END{print ind," ",sum}'
done &> ${Headeroutfile}-x.dat
for i in `seq 0 $(((${Ny}-1)))`
do 
   more ${Activityfile} | awk -v ind="$i" 'BEGIN{sum=0}{if($3==ind){sum += $1}}END{print ind," ",sum}'
done &> ${Headeroutfile}-y.dat
for i in `seq 0 $(((${Nz}-1)))`
do 
   more ${Activityfile} | awk -v ind="$i" 'BEGIN{sum=0}{if($2==ind){sum += $1}}END{print ind," ",sum}'
done &> ${Headeroutfile}-z.dat

