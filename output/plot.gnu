plot_every = 1

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.95

set xrange [-1:1]

plot_file (i, n) = sprintf ("test_angle_%d_%04d.dat", i, n)

unset key

do for [i =0:4:1] {
	set yrange [-0.1:1.1]

	set tmargin at screen 0.95
	set bmargin at screen 0.55
	
	plot for [n = 0:99:plot_every] plot_file (i, n) u 2:3 w lp
	
	set yrange [-20:20]

	set tmargin at screen 0.45
	set bmargin at screen 0.05
	
	plot for [n = 0:99:plot_every] plot_file (i, n) u 2:4 w lp
}

unset multiplot