
# CMBslow

**CMBslow** is a set of Python 3 functions for the fast calculation of the CMB power spectra as a function of the cosmological parameters $$h, \Omega_m, \text{and } \Omega_b$$. It is inspired by V. Mukhanov's impressive analytical derivation as in <https://arxiv.org/pdf/astro-ph/0303072.pdf>.

## Installation

CMBslow is not available for `pip` (yet). The main functions `IO.py` and `animate_sweep.py` however can be downloaded and plugged into any Python3 compiler directly.

## Usage

### Simple plot example
Import everything from *IO* and some utilities
```python
from IO import *
import numpy as np
import matplotlib.pyplot as plt
```
Set parameters, e.g. for the ΛCDM model
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

### Animation example
Import the animation function and needed utilities
````python
from animate_sweep import *
import numpy as np
````
Set up the conditions
````python
lmin = 200
lmax = 1000
LCDM = [0.3, 0.04,]  # parameter order [om_m, om_b, h], leave out the parameter you'd like to vary
var_range = np.linspace(0.6, 0.7, 21)
var_param = 'h'  # variation of h
````
Now call the animation function
````python
animate_sweep(minl, maxl, LCDM, var_param,  var_range)
````

# Authors

* Felix Dusel (felix.dusel@stud-mail.uni-wuerzburg.de)

This project has been developed for:
Poster Project, Cosmology ASTR407, lectured by Doulas Scott.  
University of British Columbia, Vancouver, Canada


# License
Every part of CMBslow is available under the MIT license.
