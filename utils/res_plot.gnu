set term postscript enhanced color
set out 'residue.ps'

set logscale y
set xlabel 'Number of iterations'
set ylabel 'Total residual'
p 'flo3d.res' u 1:3 w l lw 2

set logscale y
set xlabel 'Number of iterations'
set ylabel 'Residual'
p 'flo3d.res' u 1:4 t 'Density' w l lw 2, \
  'flo3d.res' u 1:5 t 'x momentum' w l lw 2, \
  'flo3d.res' u 1:6 t 'y momentum' w l lw 2, \
  'flo3d.res' u 1:7 t 'z momentum' w l lw 2, \
  'flo3d.res' u 1:8 t 'Energy' w l
