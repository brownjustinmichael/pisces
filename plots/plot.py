import IPython.html.widgets
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)
matplotlib.rc('axes', **{'linewidth': 3})

plt.style.use ("dark_background")

datasets = []
datasets.append (netCDF4.Dataset ('../sims/layers/layers_2/sc_layer_00.cdf', 'r'))
datasets.append (netCDF4.Dataset ('../sims/layers/layers_2/sc_layer_01.cdf', 'r'))
datasets.append (netCDF4.Dataset ('../sims/layers/layers_2/sc_layer_02.cdf', 'r'))
datasets.append (netCDF4.Dataset ('../sims/layers/layers_2/sc_layer_03.cdf', 'r'))

ts = {}
per = 2

xdata = []
zdata = []
tdata = []
divdata = []

for dataset in datasets:
    xdata.append (dataset.variables ["x"])
    zdata.append (dataset.variables ["z"])
    tdata.append (dataset.variables ["S"])
    divdata.append (dataset.variables ["div"])

record = -1
    
x = sum ([ma.getdata (data [record]) for data in xdata])
y = sum ([ma.getdata (data [record]) for data in zdata])
t = sum ([ma.getdata (data [record]) for data in tdata])
d = sum ([ma.getdata (data [record]) for data in divdata])

# # create a simple animation
fig = plt.figure(figsize = (6.5, 16))
ax = plt.subplot (111)
# vmin = -max (np.max (np.abs (t - y / 1.1), np.min (np.abs (t - y / 1.1))))
# vmax = max (np.max (np.abs (t - y / 1.1), np.min (np.abs (t - y / 1.1))))
p = ax.pcolormesh (x.T [::per,::per], y.T [::per,::per], t.T [::per,::per] - y.T [::per,::per] / 1.1, cmap = plt.get_cmap ("RdBu_r"))

cb = plt.colorbar (p)

ax.set_xlim ((x [0,0], x [-1,-1]))
ax.set_ylim ((y [0,0], y [-1,-1]))

ax.set_xlabel ("x", fontweight='bold', fontsize = 22)
ax.set_ylabel ("z", fontweight='bold', fontsize = 22)
cb.set_label ("Mean Molecular Weight", fontweight='bold', fontsize = 22)
# x = np.linspace(0, 10, 1000)

lines = [ax.axhline (i, -200, 200, lw = 0.25) for i in y [0,::8].flat]
vlines = [ax.axvline (i, -200, 200, lw = 0.25) for i in x [::8,0].flat]

blines = [ax.axhline (ma.getdata (i [record]) [np.logical_not (ma.getmask (i [record]))] [0], -200, 200, lw = 3.0) for i in zdata [1:]]

plt.tight_layout ()
plt.savefig ("output.pdf", transparent = True)
