import matplotlib.pyplot as plt
import numpy as np
import os

# plt.style.use ("presentation")

directory = os.path.dirname(os.path.realpath(__file__))

mw = np.genfromtxt (os.path.join (directory, "../data/moore_weiss_rb.dat"))
p = np.genfromtxt (os.path.join (directory, "../data/pisces_rb.dat"))

fig, ax = plt.subplots (1, 1)

ax.plot (mw [:,0], mw [:,1], label = "Moore & Weiss (1973)")
ax.plot (p [:,0], p [:,1], lw = 0, marker = "x", label = "ISCES (1 Element)")
ax.plot (p [:,0], p [:,2], lw = 0, marker = "o", fillstyle = "none", label = "ISCES (2 Elements)")

ax.set_xscale ("log")
ax.set_yscale ("log")

ax.legend ()

ax.set_xlabel ("$R/R_{c}$")
ax.set_ylabel ("Nu")

plt.tight_layout ()

plt.savefig ("verify.pdf")