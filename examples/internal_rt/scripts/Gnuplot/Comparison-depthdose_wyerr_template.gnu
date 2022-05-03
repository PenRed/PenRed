depthfileheaderPenR="__srcdepthPenR__"
depthfileheaderGate="__srcdepthGate__"

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
set output "Comparison-depth-dose-Z-axis-wyerr-Gate-PenRed.png"
plot 'PenRed/'.depthfileheaderPenR.'-z.dat' u 1:2:3 pointtype 7 pointsize 0.25 t ' PenRed' w yerr, 'Gate/'.depthfileheaderGate.'-z.dat' u 1:2:3 pointtype 5 pointsize 0.25 t '   Gate' w yerr
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice X"
xtop=__nx__+1
set xrange [-1:xtop]                                                                
set title "Depth dose curve along the X axis"
set output "Comparison-depth-dose-X-axis-wyerr-Gate-PenRed.png"
plot 'PenRed/'.depthfileheaderPenR.'-x.dat' u 1:2:3 pointtype 7 pointsize 0.25 t ' PenRed' w yerr, 'Gate/'.depthfileheaderGate.'-x.dat' u 1:2:3 pointtype 5 pointsize 0.25 t '   Gate' w yerr
#pause -1 'Press enter to continue'
#set logscale y
set xlabel "Slice Y"
xtop=__ny__+1
set xrange [-1:xtop]                                                                
set title "Depth dose curve along the Y axis"
set output "Comparison-depth-dose-Y-axis-wyerr-Gate-PenRed.png"
plot 'PenRed/'.depthfileheaderPenR.'-y.dat' u 1:2:3 pointtype 7 pointsize 0.25 t ' PenRed' w yerr, 'Gate/'.depthfileheaderGate.'-y.dat' u 1:2:3 pointtype 5 pointsize 0.25 t '   Gate' w yerr
#pause -1 'Press enter to continue'

