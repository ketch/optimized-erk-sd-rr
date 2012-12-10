Instructions on how to run the wedge_2d testcases

# KNOWN TO WORK WITH 
# https://github.com/wdeconinck/coolfluid3/commit/7c0d8fd89523f480189161f0cf052c7ff9ffc71b
# so download this version of coolfluid

1) prepare your shell environment:

# assuming CF3_BUILD_DIR is the directory where coolfluid is built
# assuming CF3_DEPS_DIR is the directory where third party dependencies are located

# to add a mesh-post-processing application to the path
export PATH=$CF3_BUILD_DIR/cf3/Tools/mesh_transformer:$PATH

# to find third party dependencies
export DYLD_LIBRARY_PATH=$CF3_DEPS_DIR:$DYLD_LIBRARY_PATH

# to find coolfluid in the path
export PYTHONPATH=$CF3_BUILD_DIR/dso/:$PYTHONPATH

2) generate the mesh
  - open wedge_unstr.geo with Gmsh
  - save mesh as wedge_unstr.msh in this directory

3) generate an initial condition using P0 elements (you can use as many cpu's as you like)
  - run in this directory:   
        mpirun -np 2 python wedge_2d_generate_init.py
  - consolidate distributed gmsh files into 1 file (for easy visualization and loading)
        coolfluid-mesh-transformer --input wedge_init_P*.msh --output wedge_init_consolidated.msh
        
4) run starting from initial condition (you can use as many cpu's as you like)
  - run in this directory:
        mpirun -np 2 python wedge_2d_run_from_init.py
  - This will dump every time-step of 1 second files of name wedge_t<time in seconds>.msh
    However visualizing these files with gmsh is difficult because there is 1 file for each cpu.
    Consolidate the files for better viewing, and save in any mesh-format you like.
    To visualize them for gmsh:
        coolfluid-mesh-transformer --input wedge_t1_P*.msh --output view.msh
    To visualize them for tecplot:
        coolfluid-mesh-transformer --input wedge_t1_P*.msh --output view.plt
