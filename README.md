
# CMBslow

**CMBslow** is a set of Python 3 functions for the fast calculation of the CMB power spectra as a function of the cosmological parameters $$h, \Omega_m, \text{and } \Omega_b$$. It is inspired by V. Mukhanov's impressive analytical derivation as in <https://arxiv.org/pdf/astro-ph/0303072.pdf>.

## Installation

CMBslow is not available for `pip` (yet). The main functions `IO.py` and `animate_sweep.py` however can be downloaded and plugged into any Python3 compiler directly.

## Usage

Import everything from *IO* and some utilities
```python
from IO import *
import numpy as np
import matplotlib.pyplot as plt
```
Set parameters, e.g. for the $$\Lambda \text{CDM}$$ model
```python
om_m = 0.3089
om_b = 0.0486
h = 0.6774
```
Generate a plot, for example for separated contributions
```
L = np.linspace(200, 1200, 1200-200+1)
osc, nonosc = C_l(L, om_m, om_b, h, sep_cont=True)
plt.plot(L, osc, label='Oscillating contribution')
plt.plot(L, nonosc, label='Non-oscillating contribution')
plt.plot(L, [c1 + c2 for c1, c2 in zip(osc, non)], label='result')
plt.legend()
plt.show()
```

# Authors

* Felix Dusel (felix.dusel@stud-mail.uni-wuerzburg.de)

This project has been developed for:
Poster Project, Cosmology ASTR407, lectured by Doulas Scott.  
University of British Columbia, Vancouver, Canada


# License
Every part of CMBslow is available under the MIT license.