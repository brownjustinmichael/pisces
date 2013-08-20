plot_every = 1
plot_x = 2
plot_y1 = 3
plot_y2 =4

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.95

set xrange [-1:1]

plot_file (i, n) = sprintf ("test_angle_%d_%04d.dat", i, n)

unset key

do for [i =0:3] {
	set yrange [-0.1:1.1]

	set tmargin at screen 0.95
	set bmargin at screen 0.55
	
	plot for [n = 0:99:plot_every] plot_file (i, n) u plot_x:plot_y1 w lp
	
	set yrange [-.2:.2]

	set tmargin at screen 0.45
	set bmargin at screen 0.05
	
	plot for [n = 0:99:plot_every] plot_file (i, n) u plot_x:plot_y2 w lp
}

unset multiplot