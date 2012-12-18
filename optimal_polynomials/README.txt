This directory contains instructions for generating the optimal stability polynomials
found in the paper.  However, the algorithm depends on a random number
generator (used to generate random start points for local optimization
searches), so your results will not match those in the paper exactly.

MATLAB code for generating an optimal stability function:

    p = 2;
    s = 5;
    sd.order = p;
    sd.upwindPar = 1;
    sd.KStep = 0.15;
    sd.thetaStep = 0.15;
    sd.psiStep = 0.15;
    doplot = 1;
    lam = semispectrum_2DSD_GenPatt(sd,doplot);
    x = real(lam);
    y = imag(lam);
    k = convhull(x,y);
    lam = complex(x(k),y(k));
    [h, polycoeff] = opt_poly_bisect(lam,s,p,'monomial','tol_bisect',1.e-7)


The subdirectory rk-stab-coeffs/ contains the coefficients of the optimal polynomials
used in the paper.  For instance, the file rk-stab-coeffs/sd-2.txt contains
the coefficients of all the optimal 2nd-order accurate polynomials. The first
six columns of the file contain, in order:

 - The degree of the polynomial
 - The order to which it approximates the exponential
 - The number of free coefficients in the optimization problem
 - The maximum stable step size
 - The maximum stable step size divided by the number of stages
 - The number of bisection iterations used to compute the polynomial
 
The remaining columns (from the seventh to the last) contain the polynomial
coefficients, in increasing powers of z (i.e., starting with the constant term).
