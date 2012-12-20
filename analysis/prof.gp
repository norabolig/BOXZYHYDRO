f="prof.0"
tf="txy.0"
vf="vxy.0"
dx=0.025

set terminal postscript eps enhanced font "Times" 24
set output "rhotk.eps"
set logscale xy
set key top left
set xlabel "Density (g/cc)"
set ylabel "Temperature (K)"
plot [1e-12:][30:3000] f u ($6*5.94e-7):8 lw 5 w l t"Sim", 200*(x/5.94e-11)**(1./2.4) w l t"T{/Symbol \265 r^{2.4}}"

set output "adiabaticity.eps"
unset logscale
set xlabel "R (AU)"
set ylabel "dlnT/dlnP"
plot f u 1:11 lw 5 w l t"Sim", f u 1:10 w l t"({/Symbol g}-1)/{/Symbol g}"

set terminal postscript eps enhanced color font "Times" 24
set output "motion.eps"
set xlabel "AU"
set ylabel "AU"
set title "Temperature With Gas Motion"
set size ratio 1
plot [-1:1][-1:1] tf u ($1*dx):($2*dx):3 w image t"", vf u ($1*dx):($2*dx):($3*30):($4*30) w vector lc rgb 'black' t""
