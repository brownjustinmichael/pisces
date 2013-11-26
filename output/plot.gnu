plot_every = 1

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.95

set xrange [-1:1]

plot_file (i, n) = sprintf ("output_%02i_%04d.dat", i, n)

unset key

set yrange [-0.1:10.1]

set tmargin at screen 0.95
set bmargin at screen 0.05

plot for [i = 0:1] for [n = 0:10:plot_every] plot_file (i, n) u 2:($3 + n) w lp ls n

unset multiplot