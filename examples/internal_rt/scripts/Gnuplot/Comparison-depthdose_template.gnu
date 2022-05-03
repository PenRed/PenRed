depthfiledirSet1="__srcdepthdirSet1__"
depthfiledirSet2="__srcdepthdirSet2__"
depthfileheaderSet1="__srcdepthSet1__"
depthfileheaderSet2="__srcdepthSet2__"
titleSet1="__titleSet1__"
titleSet2="__titleSet2__"

set xlabel "Slice Z"        
set ylabel "Dose (eV/g/cm**2)"                                 
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
set title "Depth dose curve along the Z axis"
set output "Comparison-depth-dose-Z-axis-".titleSet1."-".titleSet2.".png"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-z.dat' u 1:2 pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-z.dat' u 1:2 pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice X"
xtop=__nx__+1
set xrange [-1:xtop]                                                                
set title "Depth dose curve along the X axis"
set output "Comparison-depth-dose-X-axis-".titleSet1."-".titleSet2.".png"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-x.dat' u 1:2 pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-x.dat' u 1:2 pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice Y"
xtop=__ny__+1
set xrange [-1:xtop]                                                                
set title "Depth dose curve along the Y axis"
set output "Comparison-depth-dose-Y-axis-".titleSet1."-".titleSet2.".png"
plot depthfiledirSet1.'/'.depthfileheaderSet1.'-y.dat' u 1:2 pointtype 7 pointsize 0.25 t ' '.titleSet1 , depthfiledirSet2.'/'.depthfileheaderSet2.'-y.dat' u 1:2 pointtype 5 pointsize 0.25 t ' '.titleSet2
#pause -1 'Press enter to continue'

