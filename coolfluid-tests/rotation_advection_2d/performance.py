#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import sys
import math
import matplotlib
import matplotlib.ticker
matplotlib.rcParams.update({'font.size': 28})
import matplotlib.font_manager 
prop = matplotlib.font_manager.FontProperties(size=18) 
from matplotlib.ticker import MaxNLocator
import pylab
from pylab import *
pylab.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=10)
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
rc('text', usetex=True)
rc('font', family='serif')
import re

from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar


import pylab as pl
from matplotlib import ticker

def AutoLocatorInit(self):
    ticker.MaxNLocator.__init__(self, nbins=5, steps=[1, 2, 5, 10, 20])
ticker.AutoLocator.__init__ = AutoLocatorInit




def main(argv):
    # Check input arguments
    if len(argv) < 10:
        print "You must provide the order of the schemes for which the convergence plot have to generated and the name of the files that contains the error as a function of the time step. To create the legend you should also provide the name of the schemes"
        sys.exit(2)
    elif len(argv) > 10:
        print "Only the first three input arguments are used!"

    order = argv[0]

    files = tuple(argv[1:4])
    rk_names = tuple(argv[4:7])

    stages = tuple(argv[7:10])

    # Call the function in charge of creating the latex table
    do_plot(order,files,rk_names,stages)

    

def do_plot(order,file_names,rk_names,stages):
    ref_scheme_cfl, ref_scheme_eff_dt, ref_scheme_LInfinity, ref_scheme_L2, ref_scheme_cpu = read_file(order,file_names[0],stages[0])
    small_cfl, small_eff_dt, small_LInfinity, small_L2, small_cpu = read_file(order,file_names[1],stages[1])
    large_cfl, large_eff_dt, large_LInfinity, large_L2, large_cpu = read_file(order,file_names[2],stages[2])

    # Convert lists to numpy arrays
    ref_s_cfl = np.asarray(ref_scheme_cfl)
    ref_s_eff_dt = np.asarray(ref_scheme_eff_dt)
    ref_s_LInfinity = np.asarray(ref_scheme_LInfinity)
    ref_s_L2 = np.asarray(ref_scheme_L2)
    ref_s_cpu = np.asarray(ref_scheme_cpu)

    sm_cfl = np.asarray(small_cfl)
    sm_eff_dt = np.asarray(small_eff_dt)
    sm_LInfinity = np.asarray(small_LInfinity)
    sm_L2 = np.asarray(small_L2)
    sm_cpu = np.asarray(small_cpu)

    la_cfl = np.asarray(large_cfl)
    la_eff_dt = np.asarray(large_eff_dt)
    la_LInfinity = np.asarray(large_LInfinity)
    la_L2 = np.asarray(large_L2)
    la_cpu = np.asarray(large_cpu)

    ref_s_interp_cfl = np.linspace(min(ref_s_cfl),max(ref_s_cfl),6)
    z = np.polyfit(ref_s_cfl, ref_s_L2, 5)
    p = np.poly1d(z)
    ref_s_interp_L2 = p(ref_s_interp_cfl)
    z = np.polyfit(ref_s_cfl, ref_s_LInfinity, 5)
    p = np.poly1d(z)
    ref_s_interp_LInfinity = p(ref_s_interp_cfl)

    z = np.polyfit(ref_s_cfl, ref_s_cpu, 5)
    p = np.poly1d(z)
    ref_s_interp_cpu = p(ref_s_interp_cfl)

    tmp = np.linspace(max(ref_s_interp_cfl)+0.01,max(sm_cfl),6)
    sm_interp_cfl = np.concatenate((ref_s_interp_cfl,tmp))
    z = np.polyfit(small_cfl, small_L2, 5)
    p = np.poly1d(z)
    sm_interp_L2 = p(sm_interp_cfl)
    z = np.polyfit(small_cfl, small_LInfinity, 5)
    p = np.poly1d(z)
    sm_interp_LInfinity = p(sm_interp_cfl)
    z = np.polyfit(small_cfl, small_cpu, 5)
    p = np.poly1d(z)
    sm_interp_cpu = p(sm_interp_cfl)

    tmp = np.linspace(max(sm_interp_cfl)+0.01,max(la_cfl),10)
    la_interp_cfl = np.concatenate((sm_interp_cfl,tmp))
    z = np.polyfit(large_cfl, large_L2, 5)
    p = np.poly1d(z)
    la_interp_L2 = p(la_interp_cfl)
    z = np.polyfit(large_cfl, large_LInfinity, 5)
    p = np.poly1d(z)
    la_interp_LInfinity = p(la_interp_cfl)
    z = np.polyfit(large_cfl, large_cpu, 5)
    p = np.poly1d(z)
    la_interp_cpu = p(la_interp_cfl)


    # Make figure
    ############# 
    width, height = matplotlib.rcParams['figure.figsize']
    size = min(width, height)
    
    fig = plt.figure(figsize=(8.5, 6))


    host1 = fig.add_subplot(111)
    host1.plot(ref_s_cfl/float(stages[0]),ref_s_LInfinity,'ro-',ms=10,linewidth = 3, label = rk_names[0],markersize = 8)    
    host1.plot(sm_cfl/float(stages[1]),sm_LInfinity,'bs-',ms=10,linewidth = 3, label = rk_names[1], markersize = 8)
    host1.plot(la_cfl/float(stages[2]),la_LInfinity,'gv-',ms=10,linewidth = 3, label = rk_names[2], markersize = 8)
    
    #setp(host1.get_xticklabels(), visible=False)

    #host1.set_yscale('log')

    host1.set_ylabel(r'$||\varepsilon||_{L_{\infty}}$')
    host1.set_xlabel(r'$\nu/s$')
    l = plt.legend(loc = 'best',prop=prop)
    plt.grid(True)


    #host2 = fig.add_subplot(212,sharex=host1)
    #host2.plot(ref_s_cfl,ref_s_cpu,'ro-',linewidth = 3, label = rk_names[0],markersize = 8)
    #host2.plot(sm_cfl,sm_cpu,'bs-',linewidth = 3, label = rk_names[1], markersize = 8)
    #host2.plot(la_cfl,la_cpu,'gv-',linewidth = 3, label = rk_names[2], markersize = 8)


    #host2.set_xlabel(r'CFL')
    #host2.set_ylabel(r'CPU  [s]')
    #l = plt.legend(loc = 'best',prop=prop)
    #plt.grid(True)


    host1.yaxis.major.formatter.set_powerlimits((0,0))
    host1.xaxis.set_major_locator(MultipleLocator(0.01))
    host1.xaxis.set_minor_locator(MultipleLocator(0.005))
    plt.xlim(0.0,0.050000001)

    #host1.ticklabel_format(style='sci', axis='y') 
    #host2.yaxis.major.formatter.set_powerlimits((0,0)) 
    #host2.ticklabel_format(style='sci', axis='y') 


    fig.tight_layout()



    matplotlib.ticker.ScalarFormatter(useOffset=False)
    
    

    plt.draw()
    #plt.show()

    # Cut withe space
    plt.subplots_adjust(left=0.15, right=0.95, top=0.91, bottom=0.15)
    
    # Save figure
    plt.savefig("../../figures/annulus/error-cfl-selected-optimal-erk"+order+".eps")

    



def read_file(order,input_file_name,stage):
    r"""Read the data saved in file_name
    """

    if order == '2':
        w_order = '2nd-order'
    elif order == '3':
        w_order = '3rd-order'
    elif order == '4':
        w_order = '4th-order'
    elif order == '5':
        w_order = '5th-order'
    else:
        print "Order of the method unknown!"
        sys.exit(2)

    # Open file and read the content
    f_in = open('./'+w_order+'/'+input_file_name,'r');

    # Read file
    lines = f_in.readlines()

    # Close file
    f_in.close()

    # Read cfl number
    cfl_numbers = map(float, lines[1].split())

    # Compute effective time steps
    effective_dt = []
    for cfl in cfl_numbers:
        effective_dt.append(cfl/float(stage))

    # Read L_infinity norm
    LInfinity = map(float, lines[4].split())

    # Read L_2 norm
    L2 = map(float, lines[7].split())

    # Read CPU time [s]
    cpu_time = map(float, lines[10].split())
    
    # Return lists of float 
    return cfl_numbers,effective_dt,LInfinity,L2,cpu_time




if __name__ == "__main__":
  # will run only if module directly run
  main(sys.argv[1:])
else:
  # will run only if module imported
  print "I am being imported"



