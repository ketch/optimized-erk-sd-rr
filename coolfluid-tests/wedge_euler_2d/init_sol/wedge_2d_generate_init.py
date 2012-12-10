# KNOWN TO WORK WITH 
# https://github.com/wdeconinck/coolfluid3/commit/7c0d8fd89523f480189161f0cf052c7ff9ffc71b

import sys
import math
import time as systime

# Import coolfluid --> set path where to find it ( Not necessary if it is in your PYTHONPATH )
sys.path.append('path to cf3 directory/dso/')
from coolfluid import *

Core.environment().options().set('exception_outputs',False);
Core.environment().options().set('exception_backtrace',False);

simulation_order = 2         # Order of spectral difference scheme and ERK
end_time=150

############################
# Create simulation
############################
model   = Core.root().create_component('wedge_2d','cf3.solver.ModelUnsteady');
time    = model.create_time()
physics = model.create_physics('cf3.physics.NavierStokes.NavierStokes2D')
solver  = model.create_solver('cf3.sdm.SDSolver')
domain  = model.create_domain()
mesh    = domain.load_mesh(file = URI('wedge_pave.msh'), name = 'wedge2d');

### Configure physics
gamma = 1.4
R=287.05

M_inf = 0.2
p_inf = 1
rho_inf = 1
c_inf = math.sqrt(gamma*p_inf/rho_inf)
u_inf = M_inf*c_inf
T_inf = p_inf/(rho_inf*R)
rhoE_inf = p_inf/(gamma-1) + 0.5 * rho_inf * u_inf**2

Tt_inf = T_inf*(1+(gamma-1)/2*M_inf**2)
Pt_inf = p_inf*(1+(gamma-1)/2*M_inf**2)**(gamma/(gamma-1))

physics.options().set('gamma',gamma) \
                 .set('R',R)

### Configure solver
solver.options().set('mesh',mesh) \
                .set('time',time) \
                .set('solution_vars','cf3.physics.NavierStokes.Cons2D') \
                .set('solution_order',simulation_order)

### Configure timestepping
solver.access_component('TimeStepping').options().set('time_accurate',True);         # time accurate for initial stability
solver.access_component('TimeStepping').options().set('cfl','0.5');
solver.access_component('TimeStepping/IterativeSolver').options().set('nb_stages',1) # Runge Kutta number of stages

### Prepare the mesh for Spectral Difference (build faces and fields etc...)
solver.get_child('PrepareMesh').execute()

### Set the initial condition
init = solver.get_child('InitialConditions').create_initial_condition( name = 'uniform' )
functions = [ str(rho_inf), str(rho_inf*u_inf), '0.', str(rhoE_inf) ]
init.options().set('functions',functions)
init.execute();

### Create convection term
solver.get_child('DomainDiscretization').create_term(name = 'convection', type = 'cf3.sdm.navierstokes.Convection2D')
convection = solver.access_component('DomainDiscretization/Terms/convection')

### Create inlet boundary condition
bc_in   = solver.get_child('BoundaryConditions').create_boundary_condition(
   name = 'inlet', 
   type = 'cf3.sdm.navierstokes.BCSubsonicInletTtPtAlpha2D', 
   regions=[mesh.access_component('topology/inlet').uri()])
bc_in.options().set("Pt",str(Pt_inf)) \
               .set("Tt",str(Tt_inf)) \
               .set("alpha","0") \
               .set("gamma",gamma) \
               .set("R",R)
               
### Create outlet boundary condition
bc_out = solver.get_child('BoundaryConditions').create_boundary_condition(
   name = 'outlet', 
   type = 'cf3.sdm.navierstokes.BCSubsonicOutlet2D', 
   regions=[mesh.access_component('topology/outlet').uri(),
			mesh.access_component('topology/top').uri(),
			mesh.access_component('topology/bottom').uri()])
bc_out.options().set('p',p_inf)\
                .set('gamma',gamma)

### Create wall boundary condition
solver.get_child('BoundaryConditions').create_boundary_condition(
   name = 'wedge', 
   type = 'cf3.sdm.navierstokes.BCWallEuler2D', 
   regions=[mesh.access_component('topology/triangle').uri()])

#######################################
# SIMULATE
#######################################

# create post processing field
mesh.access_component('solution_space').create_field(name='post_proc',variables='U[vec],p[1],T[1],M[1],Pt[1],Tt[1],Cp[1],S[1]')
post_proc=mesh.access_component('solution_space/post_proc')
solution=mesh.access_component('solution_space/solution')

# fields to output:
fields = [
  mesh.access_component('solution_space/solution').uri(),
  mesh.access_component('solution_space/post_proc').uri(),
  mesh.access_component('solution_space/wave_speed').uri()
]

mesh.write_mesh(file=URI('wedge_load.msh'),fields=[
  mesh.access_component('solution_space/solution').uri(),
  mesh.access_component('solution_space/wave_speed').uri() ])

step=end_time
sim_time = 0
abort_simulation = False
while (sim_time<end_time and abort_simulation==False) :
	
	try :
		sim_time += step
		time.options().set('end_time',sim_time);   # limit the final time
		model.simulate()
	except(RuntimeError) :
		print "ERROR: Simulation results diverged. Aborting simulation"
		abort_simulation = True


	# POST PROCESSING

	for index in range(len(solution)):
		 rho=solution[index][0];
		 if (rho==0) :
			u=0
			v=0
			rhoE=0
			p=0
			T=0
			M=0
			Pt=0
			Tt=0
			Cp=0
			S=0
		 else :
		 	u=solution[index][1]/rho;
		 	v=solution[index][2]/rho;
		 	rhoE=solution[index][3];
		 	p=(gamma-1)*(rhoE - 0.5*rho*(u**2+v**2));
		 	T=p/(rho*R);
		 	M=math.sqrt(u**2+v**2)/math.sqrt(abs(gamma*p/rho));
		 	Pt=p+0.5*rho*(u**2+v**2);
		 	Tt=T*(1+((gamma-1)/2))*M**2;
		 	Cp=(p-p_inf)/(0.5*rho_inf*u_inf**2);
		 	S=p/(abs(rho)**gamma);

		 post_proc[index][0] = u;
		 post_proc[index][1] = v;
		 post_proc[index][2] = p;
		 post_proc[index][3] = T;
		 post_proc[index][4] = M;
		 post_proc[index][5] = Pt;
		 post_proc[index][6] = Tt;
		 post_proc[index][7] = Cp;
		 post_proc[index][8] = S;

	# OUTPUT SNAPSHOT

	mesh.write_mesh(file=URI('wedge_init'+str(sim_time)+'.msh'),fields=fields)

#####################
# FINAL OUTPUT
#####################
# Write final mesh that will be used to start up other simulations.
mesh.write_mesh(file=URI('wedge_init.msh'),fields=fields)

