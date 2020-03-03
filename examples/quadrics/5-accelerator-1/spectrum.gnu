# gnuplot script file for penEasy
#
# Last update:
#   2008-04-19 by JS

fact=1.0

filename='ImpactDetector1-th0-tallySpectrum-impdet-1.dat'

set style line 1 linetype -1 linewidth 1 pointtype 7 pointsize 0.5

set xlabel 'E   (eV)'
set ylabel ''
set xrange [*:*]
set nologscale x
set yrange [*:*]

nhist=1.0e5
energy=6.1e6/120
#energy=1.0

set terminal png size 1200,900 enhanced



set title 'Energy spectra of photons.'
set output 'plots/SpectrumPhotons.png' 
plot 'spectre.psf' u ($2*fact):5 title 'penclasses'  linestyle 2,\
    filename u 1:($6*nhist*energy) w points linestyle 1  title 'penmain'


set title 'Energy spectra of electrons.'
set output 'plots/SpectrumElectrons.png'
plot 'spectre.psf' u ($2*fact):4 title 'penclasses' linestyle 2,\
    filename u 1:($4*nhist*energy) w points linestyle 1  title 'penmain'



set title 'Energy spectra of positrons.'
set output 'plots/SpectrumPositrons.png'
plot 'spectre.psf' u ($2*fact):6 title 'penclasses' linestyle 2,\
    filename u 1:($8*nhist*energy) w points linestyle 1  title 'penmain'
  
   

set title 'Total energy spectra.'
set output 'plots/SpectrumTotal.png'
plot 'spectre.psf' u ($2*fact):3 title 'Total',\
    'spectre.psf' u ($2*fact):4 title 'Electrons',\
    'spectre.psf' u ($2*fact):5 title 'Photons',\
    'spectre.psf' u ($2*fact):6 title 'Positrons'
  


# EOF
