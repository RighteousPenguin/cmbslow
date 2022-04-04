from IO import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation as fa


def animate_sweep(l, params, sweepvar, paramrange, lmin=180, lmax=1200):

    if sweepvar == 'om_m':
        om_b, h = params
        var_latex = '$\Omega_m$'
    elif sweepvar == 'om_b':
        om_m, h = params
        var_latex = 'r$\Omega_b$'
    elif sweepvar == 'h':
        om_m, om_b = params
        var_latex = '$h$'
    else:
        print('Invalid sweep variable! Please choose either om_m, om_b, or h.')
        return None

    fig, ax = plt.subplots()
    ln, = plt.plot([], [])
    # ln2, = plt.plot([], [])  # add option sep_cont

    # plot options
    def init():
        ax.set_xlim(lmin, lmax)
        ax.set_ylim(-2, 8)
        ax.grid()
        ax.set_xlabel('$l$')
        ax.set_ylabel('$C_l$')
        return ln,

    def update(frame):
        xdata = l
        if sweepvar == 'om_m':  # poor implementation
            ydata = C_l(l, frame, om_b, h)
        elif sweepvar == 'om_b':
            ydata = C_l(l, om_m, frame, h)
        else:  # then sweepvar is h because of previous check
            ydata = C_l(l, om_m, om_b, frame)

        # ln.set_data(xdata, [ydata1 + ydata2 if not sep_cont else [ydata1, ydata2]])
        ln.set_data(xdata, ydata)  # temporary solution
        plt.title(f'{var_latex} = {np.round(frame,2)}')
        return ln,

    movie = fa(fig, update, frames=paramrange, init_func=init, blit=False)
    plt.show()
    return None


# arr = [om_m=0.3, om_b=0.04, h=0.7] just leave out the parameter you would like to sweep
LCDM = [0.04, 0.7]
var_range = np.linspace(0.1, 0.5, 20)
# simulation parameters
lmin = 160
lmax = 1400
ls = np.linspace(lmin, lmax, lmax-lmin+1)

# animate for sweeping om_m
animate_sweep(ls, LCDM, 'om_m',  var_range, lmin=lmin, lmax=lmax)
