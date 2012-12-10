# KNOWN TO WORK WITH 
# https://github.com/wdeconinck/coolfluid3/commit/77676d3bbc68a897204e9d3f3f83ba7c07f21787

import sys
import math
import time as systime

# Import coolfluid --> set path where to find it ( Not necessary if it is in your PYTHONPATH )
sys.path.append('path_to_cf3_build_dir/dso/')
from coolfluid import *

Core.environment().options().set('exception_backtrace',False)
Core.environment().options().set('exception_outputs',False)

# Import erk --> set path where to find it ( Not necessary if it is in your PYTHONPATH )
sys.path.append('../../rk-coeffs-ls')
import erk

simulation_order = 3         # Order of spectral difference scheme and ERK
rk_stages = 18               # number of stages of ERK
end_time = 200

########################################
# Runge Kutta time stepping parameters #
########################################
rk = erk.get_coeffs(simulation_order,rk_stages)
cfl = rk.cfl() # You can here also manually choose a different cfl number to become stable
print "average cfl per stage = ",rk.cfl()/rk.nb_stages

############################
# Create simulation
############################
model   = Core.root().create_component('wedge_2d','cf3.solver.ModelUnsteady');
time    = model.create_time()
physics = model.create_physics('cf3.physics.NavierStokes.NavierStokes2D')
solver  = model.create_solver('cf3.sdm.SDSolver')
domain  = model.create_domain()
mesh    = domain.load_mesh(file = URI('wedge_init_consolidated_P0.msh'), name = 'wedge2d');

## Configure physics
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
solver.access_component('TimeStepping').options().set('cfl',str(rk.cfl()));

solver.options().set('iterative_solver','cf3.sdm.RungeKuttaLowStorage3')
solver.access_component('TimeStepping/IterativeSolver').options() \
	.set('order',rk.order)\
	.set('nb_stages',rk.nb_stages)\
	.set('c',rk.c)\
	.set('delta',rk.delta)\
	.set('gamma1' ,rk.gamma1)\
	.set('gamma2' ,rk.gamma2)\
	.set('gamma3' ,rk.gamma3)\
	.set('beta',rk.beta)


### Prepare the mesh for Spectral Difference (build faces and fields etc...)
solver.get_child('PrepareMesh').execute()

mesh.print_tree()

### Set the initial condition
interpolator = model.get_child('tools').create_component('interpolator','cf3.mesh.actions.Interpolate')
interpolator.options().set('source',mesh.access_component('discontinuous_geometry/solution'))\
                      .set('target',mesh.access_component('solution_space/solution'))
interpolator.execute()

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


# Add probes in a few points.
solver.add_probe(
				name="surface",
				coordinate=[0.5,0.],
				functions=['p='+str(gamma-1)+'*(RhoE - 0.5*Rho*((RhoU[0]/Rho)^2+(RhoU[1]/Rho)^2))'],
				log_variables=['p']
				)
solver.add_probe(
				name="topshedpt",
				coordinate=[0.5,0.5],
				functions=['p='+str(gamma-1)+'*(RhoE - 0.5*Rho*((RhoU[0]/Rho)^2+(RhoU[1]/Rho)^2))'],
				log_variables=['p']
				)
solver.add_probe(
				name="x1",
				coordinate=[1.,0.],
				functions=['p='+str(gamma-1)+'*(RhoE - 0.5*Rho*((RhoU[0]/Rho)^2+(RhoU[1]/Rho)^2))'],
				log_variables=['p'],
				)
solver.add_probe(
				name="x2",
				coordinate=[2.,0.],
				functions=['p='+str(gamma-1)+'*(RhoE - 0.5*Rho*((RhoU[0]/Rho)^2+(RhoU[1]/Rho)^2))'],
				log_variables=['p'],
				)

step=1
sim_time = 0
abort_simulation = False
while (sim_time<end_time and abort_simulation==False) :
	sim_time += step
	time.options().set('end_time',sim_time);   # limit the final time
	
	try :
		model.simulate()
	except RuntimeError as detail:
		print "ERROR: Aborting simulation"
		print detail
		abort_simulation = True
		
	# POST PROCESSING

	for index in range(len(solution)):
		 rho=solution[index][0];
		 if (rho==0) :  # cannot devide by zero...
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

	mesh.write_mesh(file=URI('wedge_t'+str(sim_time)+'.msh'),fields=fields)