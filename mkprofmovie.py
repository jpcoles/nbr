import sys
#from pylab import figure, savefig, loglog, plotfile, close, show, xlim, ylim
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.pyplot import plotfile
from numpy import loadtxt, amax, abs, nan_to_num, all, exp, argmin, argwhere
from matplotlib import rc
#rc('text', usetex=True)

if len(sys.argv) > 1: files = sys.argv[1:]

Rmax = 1000

for f in files:
    print f
    
    fig = Figure()

    xs,rho,cummass = loadtxt(f, unpack=True)

    ri = argwhere(rho != 0)
    if len(ri):
        ri = ri[0]
    else:
        ri = 0

    minus3 = (xs/xs[ri])**-3 * rho[ri]
    minus2 = (xs/xs[ri])**-2 * rho[ri]
    minus1 = (xs/xs[ri])**-1 * rho[ri]

    ax = fig.add_subplot(111, title='density'); 

    fs = []
    if ri: fs += [xs, minus3, 'y-']
    if ri: fs += [xs, minus2, 'y-']
    if ri: fs += [xs, minus1, 'y-']
    w = rho != 0
    fs += [xs[w], rho[w], 'k-']
    if 1:
        ri = argmin(abs(cummass - (cummass[-1]/2)))
        rho_e = rho[ri]
        if rho_e:
            re = xs[ri]
            #print 're',re
            n = 3
            dn = 3*n - 1./3 + 0.0079/n
            ein = rho_e * exp(-dn*((xs/re)**(1./n) - 1))
            fs += [xs, ein, 'm-']

    ax.loglog(*fs)

#   ax.loglog(
#             xs, minus3, 'y-',
#             xs, minus2, 'y-',
#             xs, minus1, 'y-',
#             xs, ein,    'm-',
#             xs, cummass, 'c-',
#             #xs,rho_in,  'r-',
#             #xs,rho_out, 'b-',
#             xs,rho,     'k-'
#             );
    ax.set_xlim(0,Rmax); 
    #ax.set_ylim(10**-2, 10**10); 
    ax.set_ylim(10**-2, 10**5); 
    #ax.set_ylim(10**-20, 10**-9); 
    ax.grid(True)

    if 0:
        if 1:
            ax = fig.add_subplot(232, title='Normalized density')
            ax.loglog(
                      xs, rho/minus3, 'y-',
                      xs, rho/minus2, 'y-',
                      xs, rho/minus1, 'y-',
                      xs, rho/ein,    'm-',
                      );
            ax.set_xlim(0,1000); 
            ax.set_ylim(10**-4, 10**3);
            ax.grid(True)
        else:
            ax = fig.add_subplot(232, title='T,V,T/V')
            w = V != 0
            ax.loglog(xs,abs(T),   'r',
                      xs[w],abs(V[w]),   'b',
                      xs[w],abs(T[w]/V[w]), 'k'); 
            ax.set_xlim(0,100); 
            #ax.set_ylim(10**-10, 10**10); 
            ax.grid(True)
            #ax.legend()


    if 0:

        ax = fig.add_subplot(233, title='Radial Velocity sigma'); 
        w = sigma_r != 0
        ax.loglog(xs[w],nan_to_num(sigma_r[w]));
        ax.set_xlim(0,Rmax); 
        #ax.set_ylim(10**0, 10**3); 
        ax.grid(True)

        ax = fig.add_subplot(234, title='Tangential Velocity sigma'); 
        w = sigma_theta != 0
        ax.loglog(xs[w],sigma_theta[w]);
        ax.set_xlim(0,Rmax); 
        #ax.set_ylim(10**0, 10**3); 
        ax.grid(True)
        
        ax = fig.add_subplot(235, title='Angular Mom.'); 
        w = j!=0
        ax.loglog(xs[w],j[w]);    
        ax.set_xlim(0,Rmax); 
        #ax.set_ylim(10**-2, 10**2); 
        ax.grid(True)

        sigma3 = sigma_r*sigma_theta**2
        w = sigma3 != 0
        #X = nan_to_num(rho/(sigma_r*sigma_theta**2))
        #X = nan_to_num((sigma_r * sigma_theta**2))
        #X[X > 1e100] = 1e100

        ax = fig.add_subplot(236, title='rho/sigma^3'); 
        ax.loglog(xs[w],rho[w]/sigma3[w]);    
        ax.set_xlim(0,Rmax); 
        #ax.set_ylim(10**-7, 10**4); 
        ax.grid(True)

    canvas = FigureCanvasAgg(fig)
    canvas.print_figure("%s.png" % f)
    #show()
    #savefig("%s.png" % f)
    #close()
    
