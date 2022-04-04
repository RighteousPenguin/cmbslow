# script to plot data of PLANCK 2018

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("PowerSpect_CMB-Fit.txt")
# data = np.loadtxt("PowerSpect_CMB-TT.txt")
Ls = data[:, 0]
Dls = Ls*(Ls+1)*data[:, 1]/2/np.pi
Cls = data[:, 1]
plt.plot(Ls, Cls)
plt.xlabel('Multipole $\ell$')
plt.ylabel('$C_\ell$')
plt.grid()
plt.savefig(r"C:\Users\Felix\Desktop\posterC_l PLANCK.png", transparent=True)
plt.show()
