#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import sys
import math
import matplotlib
matplotlib.rcParams.update({'font.size': 28})
import matplotlib.font_manager 
prop = matplotlib.font_manager.FontProperties(size=18) 
from matplotlib.ticker import MaxNLocator
import pylab
pylab.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=10)
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import re


def main(argv):
    # Check input arguments
    if len(argv) < 7:
        print "You must provide the order of the schemes for which the convergence plot have to generated and the name of the files that contains the error as a function of the time step. To create the legend you should also provide the name of the schemes"
        sys.exit(2)
    elif len(argv) > 7:
        print "Only the first three input arguments are used!"

    order = argv[0]

    files = tuple(argv[1:4])
    rk_names = tuple(argv[4:7])

    # Call the function in charge of creating the latex table
    do_plot(order,files,rk_names)

    

def do_plot(order,file_names,rk_names):
    ref_scheme_h, ref_scheme_err = read_file(order,file_names[0])
    small_h, small_err = read_file(order,file_names[1])
    large_h, large_err = read_file(order,file_names[2])

    ref_curve_h = ref_scheme_h
    indices = range(0,len(ref_curve_h))
    ref_curve_error = []
    for i in indices:
        ref_curve_error.append(float(ref_curve_h[i])**float(order))
    # Make figure
    #############

    #width, height = matplotlib.rcParams['figure.figsize']
    #size = min(width, height)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = plt.gca()
    ax.plot(ref_curve_h,ref_curve_error,'r-', label = r'Reference curve $\Delta t^2$')
    ax.plot(ref_scheme_h,ref_scheme_err,'ko-', ms=10, label = rk_names[0])
    ax.plot(small_h,small_err,'-bs',ms=10, label = rk_names[1])
    ax.plot(large_h,large_err,'g^-',ms=10,label = rk_names[2])
    

    ax.set_xscale('log')
    ax.set_yscale('log')
    l = plt.legend(loc = 'upper left',prop=prop)
    plt.grid(True)
    ax.set_yticks([1e-1,1e-3,1e-5,1e-7,1e-9,1e-11])
    ax.set_xlim(1e-4, 0.25)

    ax.set_xlabel(r'$\Delta t$')
    ax.set_ylabel(r'$|\varepsilon|$')

    plt.draw()

    # Cut withe space
    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)


    plt.savefig("../../figures/convergence/convergence-selected-optimal-erk"+order+".eps")

    



def read_file(order,input_file_name):
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
    
    # Read time steps
    indices = range(1,11)

    delta_t = []
    for i in indices:
        delta_t.append(lines[i])

    # Read time steps
    indices = range(13,23)
    error = []
    for i in indices:
        error.append(lines[i])

    return delta_t,error




if __name__ == "__main__":
  # will run only if module directly run
  main(sys.argv[1:])
else:
  # will run only if module imported
  print "I am being imported"


