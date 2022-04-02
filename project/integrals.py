from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
# from main import *
#
# print(xi)
#
# # speed of sound and integral prefactors
# cs = 1/m.sqrt(3*(1+xi))
# c1 = (1-1/(3*cs**2))**2
# c2 = cs/2
# c3 = 9*cs**3/2
#
#
# def N1int(ell, x):
#     return c1*Tp(ell, x)*m.exp(-ell**2*x**2/l_f**2)/(x**2*m.sqrt(x**2-1))
#
#
# def N2int(ell, x):
#     return c2*To(ell, x)**2*m.exp(-ell**2*x**2/l_s**2)/(x**2*m.sqrt(x**2-1))
#
#
# def N3int(ell, x):
#     return c3*To(ell, x)**2*m.sqrt(x**2-1)*m.exp(-ell**2*x**2/l_s**2)/x**4

# X = np.linspace(1,)


def integrand(x, a):
    return np.cos(a*x)/np.sqrt(x**2-1)


# check formula 91
x = np.linspace(1, 100, 1000)
y = integrand(x, 1)
plt.plot(x,y)
plt.show()
[sol, err] = quad(integrand, 1, np.inf, args=(1,))
approx = np.sqrt(np.pi)*np.cos(1+np.pi/4)
print(f' sol, approx = {sol}, {approx}')




