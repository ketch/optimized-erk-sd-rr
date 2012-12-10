This directory contains all the files related to the initial solution which will
be used to restart and perform all the computations for the journal paper.


Description of the files:

- wedge_pave.geo: contains all the instruction to create the domain and the mesh
- wedge_pave.msh: contains the mesh
    
    Grid parameters:

        - Fully unstructured
        - Buffer zone for the outlet boundary
        - 16,147 quadrilaterla cells
        - min_edge_length / max_edge_length: weighted average 0.7224 
                                             min 0.1583
                                             max 0.99

- wedge_init_consolidated.msh: contains both mesh and initial solution
- wedge_init_consolidated_P0.msh: same as wedge_init_consolidated.msh




