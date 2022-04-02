from IO import *
import matplotlib.pyplot as plt

# cosmological parameters
h = 0.7
om_m = 0.3
om_b = 0.04

lmax = 1500
L = np.linspace(160, lmax, lmax-159)
rO, rN = C_l(L, sep_cont=True)

Is = np.zeros(1400-159)
ls = np.linspace(160, lmax-159, lmax-159)
# N1s = [N1(_) + N2(_) + N3(_) for _ in ls]
# A1s = [A1(_) for _ in ls]
# A2s = [A2(_) for _ in ls]
# O1s = [100/9*m.sqrt(m.pi/rho/_)*A1(_)*m.cos(_*rho+m.pi/4) for _ in ls]
# O2s = [100/9*m.sqrt(m.pi/rho/_)*A2(_)*m.cos(2*_*rho+m.pi/4)/m.sqrt(2) for _ in ls]

# for l in range(1, 101):
    # print(quad(lambda x: N1int(x, l), 1, m.inf)[0])
    # print(quad(lambda x: N1int(x, l), 1, m.inf)[1] + quad(lambda x: N2int(x, l), 1, m.inf)[1] + quad(lambda x: N3int(x, l), 1, m.inf)[1])
    # Is[l-1] = quad(lambda x: N3int(x, l), 1, m.inf)[0]  + quad(lambda x: N2int(x, l), 1, m.inf)[0] # + quad(lambda x: N3int(x, l), 1, m.inf)[0]

# plt.plot(ls, Is)
# plt.plot(ls, N1s)
# plt.xlim([180, 220])

plt.plot(L, rO, label='Oscillating contribution')
plt.plot(L, rN, label='Non-oscillating contribution')
plt.plot(L, [a + b for a, b in zip(rO, rN)], label='result')
# plt.plot(L, A1s, label='A1')
# plt.plot(L, A2s, label='A2')
# plt.plot(L, O1s, label='O1')
# plt.plot(L, O2s, label='O2')
# plt.plot(L, [a + b for a, b in zip(O1s, O2s)], label='O1+O2')
plt.xlabel('$l$')
plt.ylabel(r'$\frac{l(l+1)C_l}{(L(l+1)Cl)_{l<30}}$')
plt.grid()
plt.legend()
plt.show()

# check prefactors for O1
# fac1 = -xi*(4/(3+3*xi))**(1/4)
# fac2 = 2*m.sqrt(cs)*(1-1/(3*cs**2))
# print(f'fac1 = {fac1}, fac2 = {fac2}')
# print(f'difference = {(fac2-fac1)/fac1}')

# check prefactors for O2
# fac1 = 1/(4*m.sqrt(3+3*xi))
# fac2 = cs/2
# print(f'fac1 = {fac1}, fac2 = {fac2}')
# print(f'difference = {fac2-fac1}')
