This file contains instructions for generating the optimal Runge-Kutta methods
found in the paper.  However, the algorithm depends on a random number
generator (used to generate random start points for local optimization
searches), so your results will not match those in the paper exactly.

MATLAB code for generating an optimal method:

    p = 2;
    s = 5;
    (load polynomial coefficients from file and store them in poly_coeff)
    ind = p+1:s;
    rk = rk_opt(s,p,'erk','acc','poly_coeff_ind',ind,'poly_coeff_val',poly_coeff(ind))


Each subdirectory erk-x-N/ contain the coefficients of the optimal methods
of order N used in the paper.  For instance, the file erk-x-2/3sstar/ERK-3-2.txt contains
the coefficients of the optimal 3-stage, 2nd-order method.

The file erk.py is a convenient script for extracting the method
coefficients for manipulation in Python.
