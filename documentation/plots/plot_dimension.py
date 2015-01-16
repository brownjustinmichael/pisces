import matplotlib
import matplotlib.pyplot as plt
import numpy as np

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)
matplotlib.rc('axes', **{'linewidth': 3})

plt.style.use ("dark_background")

paddi_data = np.genfromtxt ("../data/OUT01") [::10]
pisces_data = np.genfromtxt ("../data/sc_gwave_test_PISCES", names = True)

fig, axes = plt.subplots (2, 1, sharex = True, figsize = (6.5, 6.5))

ax = axes [0]
ax.plot (paddi_data [:,1], -paddi_data [:,7], color = "white", linestyle = "--")
ax.plot (pisces_data ["t"] + 222, pisces_data ["wT"] / 10000, color = "#ff6060")

ax.set_ylabel ("Heat Flux", fontweight='bold', fontsize = 22)

ax.set_xlim ((400, 1500))
ax.set_ylim ((-2, 2))
ax.grid ()

ax = axes [1]
ax.plot (paddi_data [:,1], -paddi_data [:,8], color = "white", linestyle = "--", label = "PADDI")
ax.plot (pisces_data ["t"] + 222, pisces_data ["wS"] / 10000, color = "#6090ff", label = "PISCES")

ax.set_xlim ((400, 1500))
ax.set_ylim ((-4, 4))

ax.set_ylabel ("Comp. Flux", fontweight='bold', fontsize = 22)
ax.set_xlabel ("t (Thermal Times)", fontweight='bold', fontsize = 22)

ax.legend (loc = 'lower right')
ax.grid ()

plt.tight_layout ()
# plt.subplots_adjust (bottom = 0.15)

# plt.show ()
plt.savefig ("plot_comparison.pdf", transparent = True)