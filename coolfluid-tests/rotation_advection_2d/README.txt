Steps for simulating a rotating gaussian pulse in an annulus

The mesh
--------
- Adapt the Gmsh .geo file to the desired grid density.
- Create the mesh from this file and save it in any desired location

The simulation file
-------------------
- Adapt rotation_advection_2d.py to point to the correct coolfluid distribution
- Adapt rotation_advection_2d.py to point to the saved mesh
- Change parameters such as cfl-number, SD-order, RK-order and RK-parameters

Running
-------
- Don't forget to have the coolfluid dependencies in the DYLD_LIBRARY_PATH
- Launch using "mpirun -np 4 python rotation_advection_2d.py"   (or use different number of processes)

Results
-------
- The L2 norm of the difference between the simulated solution and exact solution is given after simulation --> only thing that really matters
- visual results are obtained with
   * tecplot: annulus_P*.plt
   * gmsh:    annulus_P*.msh

Notes
-----
Gmsh is not good for parallel visualization for now (the Gmsh people are working on it). Every file needs to be looked at separately. It is better to run with 1 cpu in case you want Gmsh output.
Tecplot can visualize parallel files together.