#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import sys
import math
import matplotlib
matplotlib.rcParams.update({'font.size': 28})
import pylab
pylab.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=10)
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import re



def main(argv):
    r"""Main function that coordinates the generation of the latex table with 
    the maximum lineraly stable time step, the relative efficiency improvment
    over a given Runge-Kutta scheme and the coefficients of the stability 
    polynomials. The first two results are printed in the same table whereas
    the coefficients are printed in a separate table because they will be added
    into the appendix of the paper.

    How to run this python script in a stand-alone mode:
    $ python table_figures_rk.py input_file_name name_reference_scheme effective_time_step_reference_scheme leading_truncation error coefficient_reference_scheme

    Example: python tables_figures_rk.py SD-5.txt 'ERKF(6,5)' 2.6916e-2 3.3557e-3 
    """

    # Check input arguments
    if len(argv) < 4:
        print "You must provide the name of the input txt file, the name of the scheme used for comparison, its maximum linearly allowable effective time step and the leading truncation error constant!"
        sys.exit(2)
    elif len(argv) > 4:
        print "Only the first four input arguments are used!"

    # Construct a tuple from the user inputs which contains the informations of
    # the reference standard explicit Runge-Kutta scheme used for comparison
    comp = tuple(argv[1:])

    # Call the function in charge of creating the latex table
    print_table(argv[0],comp)


def print_table(input_file_name,comp):
    r"""Create latex table with the coefficients of the stability polynomials 
    and the maximum lineraly stable time step.
    Inputs:
        input_file_name: name of the file containing the coefficients and the 
                         time step.

    Output: 
        .tex file containg the latex code to generate the table.
    """
    
    # Read data from input file
    stage_number,h,eff_h,coeffs,order = read_max_time(input_file_name)

    # Read leading truncation error coefficients from generated Runge-Kutta
    # coefficients
    if input_file_name.find('2')!= -1:
        path_derived_input_file_name = '../rk-coeffs-ls/erk-x2/3sstar/ERK-'
    elif input_file_name.find('3')!= -1:
        path_derived_input_file_name = '../rk-coeffs-ls/erk-x3/3sstar/ERK-'
    elif input_file_name.find('4')!= -1:
        path_derived_input_file_name = '../rk-coeffs-ls/erk-x4/3sstar/ERK-'
    elif input_file_name.find('5')!= -1:
        path_derived_input_file_name = '../rk-coeffs-ls/erk-x5/3sstar/ERK-'
    else:
        print "Order of the method unknown!"
        sys.exit(2)

    errcoeffs = read_err_coeffs(path_derived_input_file_name,stage_number,order)

    # Print latex table
    make_tables(stage_number,h,eff_h,coeffs,order,comp,errcoeffs)
    


def read_max_time(input_file_name):
    r"""Read the data saved in file_name
    """
    # Open file and read the content
    f_in = open('../rk-stab-coeffs/'+input_file_name,'r');

    # Read file
    lines = f_in.readlines()

    # Close file
    f_in.close()
    
    # Read number of stability functions
    stab_number = int(float(lines[1]))

    # Start coefficients
    start_c = 5

    h = []
    eff_h =  []
    stage_number = []
    line = lines[start_c+stab_number-1].split()
    coeffs = []

    # Read effective time step and stability polynomial coefficients 
    index = range(start_c,start_c+stab_number)
    for i in index:
        line = lines[i].split()
        stage_number.append(int(line[0]))
        h.append(float(line[3]))
        eff_h.append(float(line[4]))
        co = [float(x) for x in line[6:]]
        coeffs.append(co)
    
    # Read order of the Runge-Kutta scheme
    order = line[1]

    f_in.close()

    return stage_number,h,eff_h,coeffs,order



def read_err_coeffs(path_derived_input_file_name,stage_number,order):
    r"""Give an order of accuracy read the leding error coefficients of the 
    optimized explicit Runge-Kutta schemes
    """
    # Number of erk scheme to scan
    stab_number = len(stage_number)

    # List of errcoeffs
    errcoeffs = []

    index_erk = range(0,stab_number)
    for i in index_erk:
        # Open file, search for the leding error coefficient and read it
        derived_input_file_name = path_derived_input_file_name+order+'-'+str(stage_number[i])+'.txt'
        errcoeffs.append(find_errorcoeff(derived_input_file_name))
    
    return errcoeffs
 

def find_errorcoeff(derived_input_file_name):
    r"""Find error coefficient in the file
    """
    # Open file
    f_in = open(derived_input_file_name,'r');
    
    # Search errcoeff and read its value
    word = "errcoeff"
    line_num = 0

    file_list = f_in.readlines()
    for line in file_list:
        line_num += 1
        if line.find(word) >= 0:
            errcoeff = file_list[line_num]
            break

    # Close file
    f_in.close()

    return errcoeff
            

    
def make_tables(stage_number,h,eff_h,coeffs,order,comp,errcoeffs):
    r"""Print latex code to generate table in latex format
    """

    # Number of stability polynomial to output
    stab_number = len(stage_number)

    # Which is the maximum number of stages?
    max_s = max(stage_number)

    # Output files name 
    output_file_name_table_efficiency = 'table-efficiency-optimal-erk'+order+'.tex'

    # Open output_file_name_table_efficiency and write latex table
    f_out = open('./efficiencies/'+output_file_name_table_efficiency, 'w')

    if order == '2':
        w_order = "second"
    elif order == '3':
        w_order = 'third'
    elif order == '4':
        w_order = 'fourth'
    elif order == '5':
        w_order = 'fifth'
    elif order == '6':
        w_order = 'sixth'
    else:
        print "Order of the method unknown!"
        sys.exit(2)

    # Gain over the scheme passed as input
    gs = []
    stab_eff = []
    acc_eff = []

    # Write step size file
    ######################
    f_out.write("\\begin{table}\n")
    f_out.write("\\centering\n")
    f_out.write("\\begin{tabular}{ccc|c|c} \n")
    f_out.write("& & &\multicolumn{2}{c}{Comparison with "+comp[0]+"}\\\ \hline \n")
    f_out.write("\multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & &\\\ \n")
    f_out.write("\multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & & $\chi^{stab}$ & $\chi^{acc}$ \\\ \n")
    f_out.write("\multicolumn{1}{c|}{$s$} & \multicolumn{1}{c|}{$H/s$} & $C^{(p+1)}$ & $(H/s)_{"+comp[0]+"}$ = %.4e" % float(comp[1]))
    f_out.write(" & $(C)^{(p+1)}_{"+comp[0]+"}$ = %.4e" %float(comp[2]))
    f_out.write("\\\ \n")
    f_out.write("\multicolumn{1}{c|}{}& \multicolumn{1}{c|}{} & & & \\\ \hline \n")
    #f_out.write("& &"+comp[0]+", $H/s = $ "+comp[1]+"\\\ \hline \n")

    match = re.search(r'\d', comp[0])
    pos = match.start() if match else 'No digits found'
    number_stage_ref_rk =  float(int(comp[0][pos]))

    index_stab = range(0,stab_number)
    for i in index_stab:
        f_out.write("\multicolumn{1}{c|}{%s} &"  % stage_number[i])
        f_out.write("\multicolumn{1}{c|}{%.4e} &"  % eff_h[i])

        f_out.write(" %.4e &"  % float(errcoeffs[i]))

        # Compute relative efficiency improvment on the time step
        #perc = round(np.abs((eff_h[i]-float(comp[1]))/float(comp[1]))*100.0,2)

        gain = round(np.abs(eff_h[i]/float(comp[1])),2)
        stab_eff.append(gain)
        gs.append(stage_number[i])

        f_out.write(" %.2f &"  % stab_eff[i])

        acc_eff.append(number_stage_ref_rk/stage_number[i]*(float(comp[2])/float(errcoeffs[i]))**(1.0/(float(int(order)))))


        f_out.write(" %.2f \\\ \n"  % acc_eff[i])
        
        
    f_out.write("\hline \n")
    f_out.write("\\end{tabular} \n")
    f_out.write("\\caption{Step size, stability efficiency $\chi^{(stab)}$ and estimation of global accuracy efficiency $\chi^{(acc)}$ of the optimized "+w_order+"-order explicit Runge--Kutta methods. \n")
    f_out.write("\\label{tab:efficiencies-optimal-erk"+order+"}\n")

    f_out.write("\\end{table} \n")

    f_out.close()


    # Make figure
    ##############
    width, height = matplotlib.rcParams['figure.figsize']
    size = min(width, height)
    
    fig = plt.figure(figsize=(8.5, 6))

    # Stability efficiency
    host = fig.add_subplot(111)
    par1 = host.twinx()
    par2 = host.twinx()

    par2.spines["right"].set_position(("axes", 1.35))
    
    make_patch_spines_invisible(par2)

    par2.spines["right"].set_visible(True)

    stab, = host.plot(gs,stab_eff,'ko-',ms=10,linewidth = 3)
    acc, = par1.plot(gs,acc_eff,'bs-',ms=10, linewidth = 3)
    cfl, = par2.plot(gs,h,'r^-',ms=10, linewidth = 3)

    host.set_xlabel(r's')
    host.set_ylabel(r'$\mathbf{\chi_{stab}}$',{'color'    : 'k'})
    par1.set_ylabel(r'$\mathbf{\chi_{acc}}$',{'color'     : 'b'})
    par2.set_ylabel(r'$\nu_{stab}$',{'color'     : 'r'})

    host.yaxis.label.set_color(stab.get_color())
    par1.yaxis.label.set_color(acc.get_color())
    par2.yaxis.label.set_color(cfl.get_color())
    
    x_min = int(stage_number[0]-1)
    x_max = int(stage_number[stab_number-1]+1)
    host.set_xlim(x_min, x_max)
    host.xaxis.set_major_locator(MultipleLocator(2))
    host.xaxis.set_minor_locator(MultipleLocator(1))
    
    min_stab_eff = min(stab_eff[:])
    max_stab_eff = max(stab_eff[:])
    y_min_stab = np.around(min_stab_eff, decimals = 1) - 0.100001
    y_max_stab = np.around(max_stab_eff, decimals = 1) + 0.100001  
    host.set_ylim(1.0, 1.7000001)
    host.yaxis.set_minor_locator(MultipleLocator(0.05))

    min_stab_acc = min(acc_eff[:])
    max_stab_acc = max(acc_eff[:])
    y_min_acc = np.around(min_stab_acc, decimals = 1) - 0.1000001
    y_max_acc = np.around(max_stab_acc, decimals = 1) + 0.1000001
    par1.set_ylim(y_min_acc, y_max_acc)
    par1.yaxis.set_minor_locator(MultipleLocator(0.05))

    min_stab_cfl = min(h[:])
    max_stab_cfl = max(h[:])
    y_min_cfl = np.around(min_stab_cfl, decimals = 1) - 0.1000001
    y_max_cfl = np.around(max_stab_cfl, decimals = 1) + 0.1000001
    par2.set_ylim(y_min_cfl, y_max_cfl)
    par2.yaxis.set_minor_locator(MultipleLocator(0.1))


    plt.draw()

        
    # Cut withe space
    plt.subplots_adjust(left=0.18, right=0.68, top=0.95, bottom=0.15)

    plt.savefig("../figures/efficiencies/efficiencies-optimal-erk"+str(order)+".eps")


    # Output files name 
    output_file_name_efficiency_figure = 'figure-efficiency-optimal-erk'+order+'.tex'
        
    # Open output_file_name_figure_efficiency and write latex table
    f_out = open('../figures/efficiencies/'+output_file_name_efficiency_figure, 'w')


    f_out.write("\\begin{figure}[htbp!]\n")
    f_out.write("\\centering\n")
    
    f_out.write("\\includegraphics[scale=0.5]{./figures/efficiencies/efficiencies-optimal-erk"+order+".eps} \n")
    f_out.write("\\caption{Efficiencies of the optimized "+w_order+"-order explicit Runge--Kutta methods over the "+comp[0]+".} \n")
    f_out.write("\\label{fig:efficiencies-optimal-erk"+order+"}\n")
    f_out.write("\\end{figure}\n")



    # Write stability polynomial coefficients files
    ###############################################

    # "stab_page_wanted" stability polynomial per page
    stab_page_wanted = 6
    stab_page_needed_minus_one = int(stab_number/stab_page_wanted)

    for k in range(0,stab_page_needed_minus_one):     
        # Open output_file_name_coeffs and write latex table
        output_file_name_coeffs = 'table-stab-coeffs-optimal-erk'+order+'-'+str(k)+'.tex'
        f_out = open('./stab_poly_coeffs/'+output_file_name_coeffs, 'w')

        f_out.write("\\begin{sidewaystable}[htbp!]\n")
        f_out.write("\\centering\n")
        f_out.write("\\begin{tabular}{l|")
        f_out.write("l|"*(stab_page_wanted-1)+"l} \n")
        f_out.write(" & \multicolumn{"+str(stab_page_wanted)+"}{c}{Monomial basis coefficients} \\\ \hline \n")

        f_out.write("$s$")
        for l in range(0,stab_page_wanted):
            f_out.write(" & "+str(stage_number[l+stab_page_wanted*k]))

        f_out.write(" \\\ \hline \n ")

        index_stab = range(stab_page_wanted*k,stab_page_wanted*(k+1))
        index_coeff = range(0,stage_number[stab_page_wanted*(k+1)-1]+1)
        for j in index_coeff:
            for i in index_stab:
                if len(coeffs[i])-1<j:
                    f_out.write(" &  ")
                else:
                    if (j<int(float(order)+1)):
                        f_out.write(" & %.0e " % coeffs[i][j])
                    else:
                        f_out.write(" & %.8e " % coeffs[i][j])
                

            f_out.write(" \\\ \n ")
            

        f_out.write("\hline \n")
        f_out.write("\\end{tabular} \n")
        f_out.write("\\caption{Stability polynomial coefficients of the optimal "+w_order+" order explicit RK methods for $s="+str(stage_number[stab_page_wanted*k:stab_page_wanted*(k+1)])+"$}\n")
        f_out.write("\\label{tbl:stab-coeffs-optimal-erk"+order+"-"+str(k)+"}\n")
        f_out.write("\\end{sidewaystable} \n")
        
        f_out.close()

    
    if stab_number % stab_page_wanted != 0:

        # Remaining stability polynomials
        stab_left = stab_number-stab_page_needed_minus_one*stab_page_wanted

        # Open output_file_name_coeffs and write latex table
        output_file_name_coeffs = 'table-stab-coeffs-optimal-erk'+order+'-'+str(stab_page_needed_minus_one)+'.tex'
        f_out = open('./stab_poly_coeffs/'+output_file_name_coeffs, 'w')

        f_out.write("\\begin{sidewaystable}[htbp!]\n")
        f_out.write("\\centering\n")
        f_out.write("\\begin{tabular}{l|")
        f_out.write("l|"*(stab_left-1)+"l} \n")
        f_out.write(" & \multicolumn{"+str(stab_left)+"}{c}{Monomial basis coefficients} \\\ \hline \n")


        f_out.write("$s$")
        for l in range(0,stab_left):
            f_out.write(" & "+str(stage_number[l+stab_page_wanted*stab_page_needed_minus_one]))

        f_out.write(" \\\ \hline \n ")

        index_stab = range(stab_page_wanted*stab_page_needed_minus_one,stab_number)
        index_coeff = range(0,max_s+1)
        for j in index_coeff:
            for i in index_stab:
                if len(coeffs[i])-1<j:
                    f_out.write(" &  ")
                else:
                    if (j<int(float(order)+1)):
                        f_out.write(" & %.0e " % coeffs[i][j])
                    else:
                        f_out.write(" & %.8e " % coeffs[i][j])
            f_out.write(" \\\ \n ")
         

        f_out.write("\hline \n")
        f_out.write("\\end{tabular} \n")
        f_out.write("\\caption{Stability polynomial coefficients of the optimal "+w_order+"order RK methods for $s="+str(stage_number[stab_page_wanted*stab_page_needed_minus_one:])+"$}\n")
        f_out.write("\\label{tbl:stab-coeffs-optimal-erk"+order+"-"+str(stab_page_needed_minus_one)+"}\n")
        f_out.write("\\end{sidewaystable} \n")
        
        f_out.close()


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)
        


if __name__ == "__main__":
  # will run only if module directly run
  main(sys.argv[1:])
else:
  # will run only if module imported
  print "I am being imported"

