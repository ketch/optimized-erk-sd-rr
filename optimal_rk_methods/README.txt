This file contains instructions for generating the optimal Runge-Kutta methods
found in the paper.  However, the algorithm depends on a random number
generator (used to generate random start points for local optimization
searches), so your results will not match those in the paper exactly.

MATLAB code for generating optimal method:

    p = 2;
    s = 5;
    (load polynomial coefficients from file)
    rk = rk_opt(s,p,'erk','acc','poly_coeff_ind',???,'poly_coeff_val',???)


Each subdirectory erk-xN/ contain the coefficients of the optimal methods
of order N used in the paper.  For instance, the file erk-x2/3sstar/ERK-2-3.txt contains
the coefficients of the optimal 3-stage, 2nd-order method.

The file erk.py is a convenient script for extracting the method
coefficients for manipulation in Python.
