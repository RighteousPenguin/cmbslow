import numpy as np
import math as m


def C_l(l, om_m=0.3, om_b=0.04, h=0.7, sep_cont=False):

    # needed functions
    def get_params(xom_m, xom_b, xh):
        xxi = 17 * xom_b * xh ** 2
        # zrzeq = 7.8*10**(-2)/(om_m*h**2)  # paper
        zrzeq = 1 / 12.8 / (xom_m * xh ** 2)  # slightly more precise
        zr = xom_m ** (-0.18) / zrzeq * (m.sqrt(1 + 1 / zrzeq) - 1) ** (-2)
        xkneq = 0.72 / m.sqrt(xom_m * xh ** 2) * xom_m ** (-0.09) / 200  # *l*x
        xl_f = 1530 * m.sqrt(1 + zrzeq) * xom_m ** 0.09
        xl_s = 0.7 * xl_f * ((1 + 0.56 * xxi) / (1 + xxi) + 0.8 / (xxi * (1 + xxi)) * m.sqrt(xom_m * xh ** 2) / (
                1 + (1 + 1 / zrzeq) ** (-1 / 2)) ** 2) ** (-1 / 2)
        # rho = 0.014*om_m**0.16*h**(3.1*0.16)/(1+0.13*xi) # paper
        rho1 = 0.014 * xom_m ** 0.16 * m.sqrt(xh) / (1 + 0.13 * xxi)  # slightly more precise
        rho2 = xom_m ** (-0.09) / m.sqrt(3 * zr * xxi) * m.log(
            (m.sqrt(xxi * (1 + zrzeq)) + m.sqrt(1 + xxi)) / (1 + m.sqrt(xxi * zrzeq)))  # even more precise
        xrho = rho1

        # speed of sound and integral prefactors
        xcs = 1 / m.sqrt(3 * (1 + xxi))
        xc1 = (1 - 1 / (3 * xcs ** 2)) ** 2
        xc2 = xcs / 2
        xc3 = 9 * xcs ** 3 / 2
        return xxi, xkneq, xl_f, xl_s, xrho, xcs, xc1, xc2, xc3

    [xi, kneq, l_f, l_s, rho, cs, c1, c2, c3] = get_params(om_m, om_b, h)

    def P(xell, xom_m, xh):
        return m.log((xom_m ** (-0.09) * xell / 200) / m.sqrt(xom_m * xh ** 2))

    def Tp(xell, x):
        if xell > 199:
            return 0.74 - 0.25*(P(xell, om_m, h) + m.log(x))
        else:
            return 1 / 4 * m.log(14 / (kneq * x * xell))

    def To(xell, x):
        if xell > 199:
            return 0.5 + 0.36*(P(xell, om_m, h) + m.log(x))
        else:
            return 0.36 * m.log(5.6 * kneq * x * xell)

    def A1(xell):
        # return 0.1*xi*((P(ell)-0.78)**2-0.43)/(1+xi)**(1/4)*m.exp(ell**2/2*(1/l_s**2 - 1/l_f**2))  # less precise
        # return 0.1*xi*((P(ell)-0.78)**2-0.43)/(1+xi)**(1/4)*m.exp(-ell**2/2*(1/l_s**2 + 1/l_f**2))  # less precise, mine
        # return -(4/(3*(1+xi)))**(1/4)*xi*Tp(ell, 1)*To(ell, 1)*m.exp(ell**2/2*(1/l_s**2 - 1/l_f**2))  # precise, paper
        return -(16 / (3 * (1 + xi))) ** (1 / 4) * xi * Tp(xell, 1) * To(xell, 1) * m.exp(
            -xell ** 2 / 2 * (1 / l_s ** 2 + 1 / l_f ** 2))  # precise, mine

    def A2(xell):
        # return 0.14*(0.5+0.36*P(ell))**2/(1+xi)**0.5  # less precise
        # return 0.14*(0.5+0.36*P(ell))**2/(1+xi)**0.5*m.exp(-ell**2/l_s**2)  # less precise, mine
        # return To(ell, 1)**2/(4*m.sqrt(3*(1+xi)))  # precise, paper
        return To(xell, 1) ** 2 / (2 * m.sqrt(3 * (1 + xi))) * m.exp(-xell ** 2 / l_s ** 2)  # precise, mine

    def O(xell):
        # return m.sqrt(m.pi/rho/ell)*(A1(ell)*m.cos(rho*ell+m.pi/4)+A2(ell)*m.cos(2*ell*rho+m.pi/4)) # paper
        return m.sqrt(m.pi / rho / xell) * (
                A1(xell) * m.cos(rho * xell + m.pi / 4) + A2(xell) * m.cos(2 * xell * rho + m.pi / 4) / m.sqrt(2))  # mine

    def N1(xell, xom_m, xh):
        # re = np.zeros(ell)
        # for indl, ll in enumerate(ell):
        #     0.063 * xi ** 2 * (P(ll) - 0.22 * (ll / l_f) ** 0.3 - 2.6) ** 2 / (1 + 0.65 * (ll / l_f) ** 1.4) * m.exp(
        #         -ll ** 2 / l_f ** 2)
        # return re
        return 0.063 * xi ** 2 * (P(xell, xom_m, xh) - 0.22 * (xell / l_f) ** 0.3 - 2.6) ** 2 / (
                1 + 0.65 * (xell / l_f) ** 1.4) * m.exp(-xell ** 2 / l_f ** 2)

    def N2(xell, xom_m, xh):
        return 0.037 / m.sqrt(1 + xi) * (P(xell, xom_m, xh) - 0.22 * (xell / l_s) ** 0.3 + 1.7) ** 2 / (
                1 + 0.65 * (xell / l_s) ** 1.4) * m.exp(-xell ** 2 / l_s ** 2)

    def N3(xell, xom_m, xh):
        return 0.033 * (1 + xi) ** (-3 / 2) * (P(xell, xom_m, xh) - 0.5 * (xell / l_s) ** 0.55 + 2.2) ** 2 / (
                1 + 2 * (xell / l_s) ** 2) * m.exp(-xell ** 2 / l_s ** 2)

    def N1int(x, xell):
        # re = np.zeros(shape=(len(x), len(ell)))
        # for indx, xx in enumerate(x):
        #     for indl, ll in enumerate(ell):
        #         re[indx, indl] = c1*Tp(ll, xx)*m.exp(-ll**2*xx**2/l_f**2)/(xx**2*m.sqrt(xx**2-1))
        # return re
        return c1 * Tp(xell, x) * m.exp(-xell ** 2 * x ** 2 / l_f ** 2) / (x ** 2 * m.sqrt(x ** 2 - 1))

    def N2int(x, xell):
        return c2 * To(xell, x) ** 2 * m.exp(-xell ** 2 * x ** 2 / l_s ** 2) / (x ** 2 * m.sqrt(x ** 2 - 1))

    def N3int(x, xell):
        return c3 * To(xell, x) ** 2 * m.sqrt(x ** 2 - 1) * m.exp(-xell ** 2 * x ** 2 / l_s ** 2) / x ** 4

    resultO = np.zeros_like(l)
    resultN = np.zeros_like(l)
    for n, ell in enumerate(l):
        resultO[n] = 100/9*(O(ell))
        resultN[n] = 100/9*(N1(ell, om_m, h) + N2(ell, om_m, h) + N3(ell, om_m, h))

    if not sep_cont:
        return resultO + resultN
    else:  # return contributions separately
        return resultO, resultN

