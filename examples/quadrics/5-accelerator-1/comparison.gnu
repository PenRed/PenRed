# gnuplot script file for penEasy
#
# Last update:
#   2008-04-19 by JS

# Factor from eV to keV:
f = 1.0e-3
sigma=0.666667

set style line 1 linetype -1 linewidth 1 pointtype 7 pointsize 0.5

set xlabel 'E   (eV)'
set ylabel 'Prob (1/eV*particle)'
set xrange [*:*]
set nologscale x
set yrange [*:*]

set title 'Energy spectra of incident particles.'

set terminal png size 1200,900 enhanced




set output 'plots/SpectrumComparisonElectrons.png'
plot 'spc-impdet-01.dat' u 1:4:($5*sigma) w yerrorbars linestyle 1 title 'penmain',\
     "ImpactDetector1-th0-tallySpectrum-impdet-1.dat" using 1:4:5 w yerrorbars linestyle 2 title 'penclasses'
     
     
set output 'plots/SpectrumComparisonPhotons.png'

plot 'spc-impdet-01.dat' u 1:6:($7*sigma) w yerrorbars linestyle 1 title 'penmain',\
     "ImpactDetector1-th0-tallySpectrum-impdet-1.dat" using 1:6:7 w yerrorbars linestyle 2 title 'penclasses'
     
     
set output 'plots/SpectrumComparisonPositrons.png'

plot 'spc-impdet-01.dat' u 1:8:($9*sigma) w yerrorbars linestyle 1 title 'penmain',\
     "ImpactDetector1-th0-tallySpectrum-impdet-1.dat" using 1:8:9 w yerrorbars linestyle 2 title 'penclasses'
     
set output 'plots/SpectrumComparisonTotal.png'

plot 'spc-impdet-01.dat' u 1:2:($3*sigma) w yerrorbars linestyle 1 title 'penmain',\
     "ImpactDetector1-th0-tallySpectrum-impdet-1.dat" using 1:2:3 w yerrorbars linestyle 2 title 'penclasses'
     
     

# EOF
