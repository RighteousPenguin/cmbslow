import camb
import numpy as np
import matplotlib.pyplot as plt
from IO import C_l


lmax = 1200
# my script to plot data of PLANCK 2018
data = np.loadtxt("PowerSpect_CMB-Fit.txt")
# data = np.loadtxt("PowerSpect_CMB-TT.txt")
Ls = data[0:lmax, 0]
Dls = Ls*(Ls+1)*data[0:lmax, 1]/2/np.pi
Cls = data[0:lmax, 1]
# plt.xlabel('Multipole $l$')
# plt.ylabel('$C_l$')
# plt.grid()
# plt.savefig(r"C:\Users\Felix\Desktop\posterC_l PLANCK.png",
#             transparent=True)
# plt.show()

# Set up a new set of parameters for CAMB
pars = camb.CAMBparams()

# This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
om_b = 0.04
h = 0.675
H0 = 100*h
pars.set_cosmology(H0=H0)  # H0=H0, ombh2=om_b*h**2, omch2=0, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=1, r=0)
pars.set_for_lmax(lmax-50, lens_potential_accuracy=0)

# calculate results for these parameters
results = camb.get_results(pars)

# get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')

# plot the total lensed CMB power spectra versus unlensed, and fractional difference
totCL = powers['total']
unlensedCL=powers['unlensed_scalar']

# Python CL arrays are all zero based (starting at L=0), Note L=0,1 entries will be zero by default.
# The different CL are always in the order TT, EE, BB, TE (with BB=0 for unlensed scalar results).
ls = np.arange(totCL.shape[0])
plt.plot(ls[199:], unlensedCL[199:, 0], color='r', label=r'CAMB for $\Lambda$CDM, unlensed')
plt.plot(Ls[199:], Cls[199:], label='PLANCK 2018')
plt.plot(Ls[199:], C_l(Ls, h=h)[199:]*1000, label='CMBslow')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_\ell$')
plt.grid()
plt.legend()
plt.savefig(r"C:\Users\Felix\Desktop\cmbslow\figures\comparison.png", transparent=True)
plt.show()
