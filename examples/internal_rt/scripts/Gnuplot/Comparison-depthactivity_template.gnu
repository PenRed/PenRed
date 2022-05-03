depthfiledirSet1="__srcdepthdirSet1__"
depthfiledirSet2="__srcdepthdirSet2__"
depthfileheaderSet1="__srcdepthSet1__"
depthfileheaderSet2="__srcdepthSet2__"
nhisSet1=__nhisSet1__
nhisSet2=__nhisSet2__
titleSet1="__titleSet1__"
titleSet2="__titleSet2__"
set xlabel "Slice Z"        
set ylabel "Probability (%)"                                 
#set style data linespoints
#set key box
#set key fixed right top vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid
set key right top vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid
#set logscale y
#set yrange [-0.005,1]
xtop=__nz__+1
set xrange [-1:xtop]                                                                
#unset logscale  
#set terminal pngcairo size 640,480
#Options are ' background "#ffffff" enhanced fontscale 1.0 size 640, 480 '
set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 640, 480
set output "Comparison-depth-activity-Z-axis-".titleSet1."-".titleSet2.".png"
set title "Depth activity curve along the Z axis"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-z.dat' u 1:($2/nhisSet1*100) pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-z.dat' u 1:($2/nhisSet2*100) pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice X"
xtop=__nx__+1
set xrange [-1:xtop]
set title "Depth activity curve along the X axis"
set output "Comparison-depth-activity-X-axis-".titleSet1."-".titleSet2.".png"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-x.dat' u 1:($2/nhisSet1*100) pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-x.dat' u 1:($2/nhisSet2*100) pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice Y"
xtop=__ny__+1
set xrange [-1:xtop]
set title "Depth activity curve along the Y axis"
set output "Comparison-depth-activity-Y-axis-".titleSet1."-".titleSet2.".png"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-y.dat' u 1:($2/nhisSet1*100) pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-y.dat' u 1:($2/nhisSet2*100) pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'

