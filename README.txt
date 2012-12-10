ABOUT THIS REPOSITORY
================

This is the reproducibility repository for "Optimized explicit Runge-Kutta schemes for the spectral difference method applied to wave propagation problems", 
a manuscript submitted in 2012 by Parsani et al. to the SIAM Journal for
Scientific Computing (SISC).  This repository serves two purposes:

- Establish and maintain provenance for all data, charts, and tables presented in the manuscript.  

In particular, the raw data for all visualizations in the paper and their
source is archived here along with, where possible, the scripts used to
generate their manuscript-ready forms.  Preference is given
to minimizing the size of the repository while providing complete data sets and
provenance for all data visualizations in the manuscript.  Consequently,
compact post-processed data sets are the preferred archival format for this
repository.  For example, if a simulation contains many time steps before
arriving at a final solution, and computational gauges are only monitoring a
small slice of the simulated domain, it is only expected that the raw data
actually used directly used in generating the figures themselves is archived.

- Assist and automate reproducibility.

The hardware/software environment along with the software used to produce every
experimental result reported in the manuscript is described in this repository.
On systems with the prerequisite software installed, the Makefile in the root
directory documents and automates reproduction of the experiments described in
the paper.  

Automation of this repository relies on MATLAB (including the Optimization and 
Global Optimization toolboxes), Python (including NumPy), as well as the packages
RK-opt, Nodepy, and Coolfluid 3, which can be found on Github.

Finally, the REPRODUCED file maintains a list of readers who were able to
reproduce the experiments in this repository, with notes on any
discrepancies found.  

If you reproduce any of this data, please let me or another author of this
paper know so we can add your results to the REPRODUCED file. 

David Ketcheson, King Abdullah University of Science & Technology
david.ketcheson@kaust.edu.sa

