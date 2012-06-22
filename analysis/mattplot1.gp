f="dxy.30mearth"
g="vxy.30mearth"

set terminal postscript enhanced color font "Helvetica" 22
set output 'adiabatic_faceon.eps'  

r0=5.2*(30e-3/320./3.)**(1./3.)
set object 1 circle at 0,0 size r0 front fs empty
set object 2 circle at 0,0 size r0/2 front fs empty
set size square
set logscale cb
set cbrange [6e-13:6e-8]
set xlabel 'AU'
set ylabel 'AU'
set title "30 M_{/Symbol \305} Planet\n Log Density (cgs) with Gas Velocity Vectors"

plot [-0.3:0.3][-0.3:0.3] f u 1:2:($3*5.94e-7) w image, g u 1:2:($3*0.5):($4*0.5) w vector lc rgb 'black' t''

