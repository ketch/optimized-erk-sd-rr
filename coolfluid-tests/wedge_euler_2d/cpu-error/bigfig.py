import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import pylab
import matplotlib.ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')



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

    # Read CPU time [s]
    cpu_time = map(float, lines[7].split())
    
    # Return lists of float 
    return cfl_numbers,effective_dt,LInfinity,cpu_time

def main(argv):
    matplotlib.rcParams.update({'font.size': 25})

    cpu = []
    err = []

    orders = ['2','2','2','3','3','3','4','4','4','5','5','5']
    files = ['wedge_2d_mid_point.txt','wedge_2d_optimal_erk-3-2.txt','wedge_2d_optimal_erk-8-2.txt','wedge_2d_heun3.txt','wedge_2d_optimal_erk-5-3.txt','wedge_2d_optimal_erk-17-3.txt','wedge_2d_kutta44.txt','wedge_2d_optimal_erk-9-4.txt','wedge_2d_optimal_erk-18-4.txt','wedge_2d_fehlberg.txt','wedge_2d_optimal_erk-10-5.txt','wedge_2d_optimal_erk-20-5.txt']
    stages = [2,3,8,3,5,17,4,9,18,6,10,20]
    labels=[str(stage)+','+str(order) for stage,order in zip(stages,orders)]
    styles = ['wo','ko','ko','wo','ko','ko','wo','ko','ko','wo','ko','ko','wo','ko','ko']

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for order, filename, stage, style in zip(orders,files, stages, styles):
        cfl, dt, this_err, this_cpu = read_file(order,filename,stage)
        cpu.append(this_cpu[-1]); err.append(this_err[-1])
        ax.semilogy(this_cpu[-1],this_err[-1],style,markersize=10)
        plt.hold(True)

    for label, x, y in zip(labels, cpu, err):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (5, 10),
            textcoords = 'offset points', ha = 'left', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),fontsize=18)
            #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

    plt.xlabel('CPU [s]',fontsize=25)
    plt.ylabel(r'$||\varepsilon||_{L_{\infty}}$',fontsize=25)
    plt.xlim(35,72)
    plt.draw()
    plt.subplots_adjust(left=0.18, right=0.95, top=0.91, bottom=0.12)
    
    plt.savefig("../../../figures/wedge/error-cpu.pdf")

    cpu_ref = []
    cpu_small = []
    cpu_large = []
    lab_ref = []
    lab_small = []
    lab_large = []
    cpu_large = []

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ind = np.arange(4)*2
    width = 0.4
    cpu_ref.append(cpu[0])
    cpu_ref.append(cpu[3])
    cpu_ref.append(cpu[6])
    cpu_ref.append(cpu[9])
    rects_ref = ax1.bar(ind+0.5, cpu_ref, width, color='w',hatch='\\')
    plt.hold(True)

    lab_ref.append(labels[0])
    lab_ref.append(labels[3])
    lab_ref.append(labels[6])
    lab_ref.append(labels[9])

    for label, x, y in zip(lab_ref,ind+0.55, cpu_ref):
        ax1.annotate(
        label, 
        xy = (x, y), xytext = (5, 10),
        textcoords = 'offset points', ha = 'left', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),fontsize=18)

    cpu_small.append(cpu[1])
    cpu_small.append(cpu[4])
    cpu_small.append(cpu[7])
    cpu_small.append(cpu[10])
    rects_small = ax1.bar(ind+width+0.5, cpu_small, width, color='w')

    lab_small.append(labels[1])
    lab_small.append(labels[4])
    lab_small.append(labels[7])
    lab_small.append(labels[10])

    for label, x, y in zip(lab_small,ind+width+0.55, cpu_small):
        ax1.annotate(
        label, 
        xy = (x, y), xytext = (5, 10),
        textcoords = 'offset points', ha = 'left', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),fontsize=18)

    cpu_large.append(cpu[2])
    cpu_large.append(cpu[5])
    cpu_large.append(cpu[8])
    cpu_large.append(cpu[11])
    rects_small = ax1.bar(ind+width*2+0.5, cpu_large, width, color='w',hatch='//')

    lab_large.append(labels[2])
    lab_large.append(labels[5])
    lab_large.append(labels[8])
    lab_large.append(labels[11])

    for label, x, y in zip(lab_large[1:4],ind[1:4]+width*2+0.55, cpu_large[1:4]):
        ax1.annotate(
        label, 
        xy = (x, y), xytext = (5, 10),
        textcoords = 'offset points', ha = 'left', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),fontsize=18)

    ax1.annotate(
    lab_large[0], 
    xy = (ind[0]+width*2.4+0.55, cpu_large[0]), xytext = (5, 10),
    textcoords = 'offset points', ha = 'left', va = 'bottom',
    bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),fontsize=18)

    #ax1.axes.get_xaxis().set_ticks([])
    plt.ylabel('CPU [s]',fontsize=25)
    plt.xlim((0,8.35))    
    plt.ylim((0,75))
    ax1.yaxis.set_major_locator(MultipleLocator(10))
    ax1.yaxis.set_minor_locator(MultipleLocator(1))

    shift = 0.70
    pylab.xticks([ind[0]+width+shift,ind[1]+width+shift,ind[2]+width+shift,ind[3]+width+shift], [r'2$^\mathrm{nd}$-order',r'3$^\mathrm{rd}$-order',r'4$^\mathrm{th}$-order',r'5$^\mathrm{th}$-order'],fontsize=23)

    ax1.xaxis.set_ticks_position('none')

    plt.draw()    

    plt.subplots_adjust(left=0.11, right=0.95, top=0.91, bottom=0.12)
    
    plt.savefig("../../../figures/wedge/cpu.pdf")



if __name__ == "__main__":
  main(sys.argv[1:])
else:
  print "I am being imported"

