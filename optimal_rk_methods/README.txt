This file contains instructions for generating the optimal Runge-Kutta methods
found in the paper.  However, the algorithm depends on a random number
generator (used to generate random start points for local optimization
searches), so your results will not match those in the paper exactly.

Here is sample MATLAB code for generating an optimal method with order 2 and 5 stages:

    p = 2;
    s = 5;
    % Load the stability polynomial coefficients 
    poly_coeff = load_poly('../optimal_polynomials/polynomial-coefficients/sd-2.txt',5)
    rk = rk_opt(s,p,'erk','acc','poly_coeff_ind',p+1:s,'poly_coeff_val',poly_coeff)


Alternatively, you can manually copy and paste the polynomial coefficients, as follows:

    p = 2;
    s = 5;
    % The following array contains the optimal stability polynomial stored in
    % optimal_polynomials/polynomial-coefficients/sd-2.txt
    poly_coeff = [1.0000000000000000E+00,	1.0000000000000000E+00,	5.0000000000000000E-01,	1.3876772158093043E-01,	2.1520726214707169E-02,	1.4956929981860788E-03];
    ind = p+1:s;
    rk = rk_opt(s,p,'erk','acc','poly_coeff_ind',ind,'poly_coeff_val',poly_coeff(ind+1))


Each subdirectory erk-x-N/ contain the coefficients of the optimal methods
of order N used in the paper.  For instance, the file erk-x-2/3sstar/ERK-5-2.txt contains
the coefficients of the optimal 5-stage, 2nd-order method.

The file erk.py is a convenient script for extracting the method
coefficients for manipulation in Python.
