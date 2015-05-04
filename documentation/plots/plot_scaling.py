import numpy as np
import matplotlib.pyplot as plt

plt.style.use ("presentation")

data = np.genfromtxt ("../data/scaling", names = True)

s = plt.scatter (data ["nmpi"], data [0] ["med"] / (data ["med"] * (data ["nmpi"] * data ["nmp"])), c = data ["nmp"], s = data ["ntt"] * 50, alpha = 0.75)

# plt.xscale ("log")
plt.xlabel ("Number of MPI Processes")
# plt.yscale ("log")
plt.ylim ([0., None])
plt.ylabel ("Strong Scaling Efficiency")
cb = plt.colorbar (s)
cb.set_ticks ([1, 4, 8])
cb.set_label ("Number of MP Threads")

plt.savefig ("plot_scaling.pdf")
