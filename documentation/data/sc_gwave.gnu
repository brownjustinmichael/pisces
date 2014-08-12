set terminal pdf enhanced color font "Times,40" size 10,5 dashed
set border linewidth 4
set object 1 rectangle from screen 0,0 to screen 1,1 behind fs solid noborder

set output "sc_gwave.pdf"

set xlabel "t (thermal diffusion times)"
set ylabel "Turbulent Flux"

set xrange [400:1500]
set yrange [-4:2]

set key out right

plot "OUT01" ev 10 u 2:(-$8) w l ls 1 lc 1 lw 4 title "PADDI, wT", \
"sc_gwave_test_PISCES" u ($1+222):($3/10000) w l ls 1 lc 2 lw 4 title "PISCES, wT", \
"OUT01" ev 10 u 2:(-$9 - 2) w l ls 3 lc 1 lw 4 title "PADDI, wS - 2", \
"sc_gwave_test_PISCES" u ($1+222):($4/10000-2) w l ls 3 lc 2 lw 4 title "PISCES, wS - 2"