set term postscript enhanced
set out 'residue.ps'

set logscale y
set xlabel 'Number of iterations'
p 'residue.dat' u 1:3 t 'Total' w l

set logscale y
set xlabel 'Number of iterations'
p 'residue.dat' u 1:4 t 'Density' w l, \
  'residue.dat' u 1:5 t 'x momentum' w l, \
  'residue.dat' u 1:6 t 'y momentum' w l, \
  'residue.dat' u 1:7 t 'z momentum' w l, \
  'residue.dat' u 1:8 t 'Energy' w l
