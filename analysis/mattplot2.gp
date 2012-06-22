f="dyz.30mearth"

set terminal postscript enhanced color font "Helvetica" 22
set output 'adiabatic_edge.eps'  

set size square
set logscale cb
set cbrange [6e-13:6e-8]
set xlabel 'AU' 
set ylabel 'AU'
set xtics out rotate by -45
set title "30 M_{/Symbol \305} Planet\n Log Density (cgs)"

plot [-0.075:0.075][-0.075:0.075] f u 1:2:($3*5.94e-7) w image t''

