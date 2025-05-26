------------------------------------------------------------------------------------------
-- Navier-Stokes equation, 2d
-- Discretization: Vertex-centered, stabilized
------------------------------------------------------------------------------------------

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

-- Required Plugins
PluginRequired("NavierStokes")
PluginRequired("ConvectionDiffusion")


------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Physical parameters

nu_a 	= util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity")
rho_a 	= util.GetParamNumber("-rho_a", 1.2, "Air Density")
inflow		= util.GetParamNumber("-inflow", 1.0, "max. inflow velocity")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")
c_init		= util.GetParamNumber("-initial concentration", 1.0, "max volume fraction")

-- Numerical parameters of the discretization
numRefs 	= util.GetParamNumber("-numRefs", 4, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 2, "number of prerefinements (parallel)")
bNoLaplace 	= util.HasParamOption("-noLaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactJac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type (no, full, weighted, lps, pos, reg)")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")

dt = util.GetParamNumber("-dt",10, "TimeStep")
timeMethod = util.GetParam("-timeMethod","euler")
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 10	)
outputFactor     = util.GetParam("-output", 1, "output every ... steps")


-- Parameters of the solver
ilu_beta	= util.GetParamNumber("-iluBeta", -0.2, "choose a negative value depending on the convection rate")


------------------------------------------------------------------------------------------
-- Geometry data and export file name
------------------------------------------------------------------------------------------

dim 		= 2
gridName = "grids/channel_40x15.ugx"
file_name ="Transport"

-- Subsets used in the problem
allSubsets = "Inner,Left, Right,Top, Bottom,"


------------------------------------------------------------------------------------------
-- Print the parameters
------------------------------------------------------------------------------------------

print (" Geometry: " .. gridName .. ", dim = " .. dim)
print (" Physical parameter:")
print ("    visc		= " .. nu_a)
print ("    inflow		= " .. inflow)
print ("    Stokes		= " .. tostring (bStokes))
print (" Numerical parameter of the discretization:")
print ("    numRefs		= " .. numRefs)
print ("    numPreRefs	= " .. numPreRefs)
print ("	noLaplace	= " .. tostring (bNoLaplace))
print ("	exactJac	= " .. tostring (bExactJac))
print ("    PecletBlend = " .. tostring (bPecletBlend))
print ("    upwind		= " .. upwind)
print ("    pac			= " .. tostring (bPac))
print ("    stab		= " .. stab)
print ("    difflength	= " .. diffLength)
print (" Numerical parameter of the solvers:")
print ("    iluBeta		= " .. ilu_beta)



------------------------------------------------------------------------------------------
-- Initialize UG4, load, refine and distribute the grid
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", dim +2))

if dim == 3 then
	fct_cmp_tbl = {"u", "v", "w", "p"}
	vel_cmp_tbl = {"u", "v", "w"}
else
	fct_cmp_tbl = {"u", "v", "p"}
	vel_cmp_tbl = {"u", "v"}
end	

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs)
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)

approxSpace:add_fct("u", "Lagrange",1,allSubsets)
approxSpace:add_fct("v", "Lagrange",1,allSubsets)
approxSpace:add_fct("p", "Lagrange",1,allSubsets)
approxSpace:add_fct("c", "Lagrange",1,allSubsets)


approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace


------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------


-- Order the DoFs:
OrderLex (approxSpace,  "x")
--OrderCuthillMcKee(approxSpace,true)
-- grid function for the solution
u = GridFunction(approxSpace)
u:set(0)

--------------------------------
--------------------------------
-- Lua Functions
--------------------------------
--------------------------------

	
---------------------------------------------------------------------- Initial Velocity

function ValueX(x,y) 
	hh=15
	return inflow*( (hh - y) * (4*y ) / (hh * hh))
end


function ValueY(x,y) 
	return 0.0
end
function ValueZ(x,y) 
	return 0
end
---------------------------------------------------------------------- Initial Pressure

function ValueP(x,y) 
	return 0.0
end
-------------------------------------------------------------- Initial VolumeFraction 

function ValueC(x,y)
	return 0.0
end

---------------------------------------------------------------------- Boundary Condition  

----------------------------------------------------------- Inlet
function BoundaryVolumeFraction(x,y)
	return c_init
end


function inflowVel2d(x, y, t)

	return ValueX(x,y),ValueY(x,y)
end


------------------------------------------------------------------------------------------
-- Compose the discretization of NavierStokes Equation
------------------------------------------------------------------------------------------

-- inner space
NavierStokesDisc = NavierStokesFV1(fct_cmp_tbl, {"Inner"})

NavierStokesDisc:set_peclet_blend(bPecletBlend) 
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_laplace(bNoLaplace)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)


----------------------------------------------------------------------- Set  inputs

NavierStokesDisc:set_density(rho_a)
NavierStokesDisc:set_kinematic_viscosity(nu_a)



------------------------------------------------- Set  NS boundary conditions

-- boundary condition at the inlet
InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inflowVel2d", "Left")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Right")

-- boundary condition at the  walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add (" Bottom,Top")

print("Navier Stokes Equation created.")

------------------------------------------------------------------------------------------
-- Compose the discretization of Transport Equation
------------------------------------------------------------------------------------------

TransportEq = ConvectionDiffusion("c", "Inner", "fv1")
TransportEq:set_velocity(NavierStokesDisc:velocity_ip())
TransportEq:set_upwind(UpwindFV1(upwind)) 


------------------------------------------------- Set   boundary conditions


-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add(BoundaryVolumeFraction, "c", "Left")
dirichletBND:add(0.0, "c", "Top,Bottom")


OutflowBND = ConvectionDiffusionOutflowFV1(TransportEq)
OutflowBND:add( "Right")	


print("Transport Equation created.")


---------------------------------------------------------------------------------------
-- Domain Disc
---------------------------------------------------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(WallDisc)
domainDisc:add(InletDisc)
domainDisc:add(OutletDisc)


domainDisc:add(TransportEq)
domainDisc:add(OutflowBND)
domainDisc:add(dirichletBND)

-- create operator from discretization


-- create time discretization
if timeMethod=="cn" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.5) -- Crank-Nicolson method
end
if timeMethod=="euler" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1) -- implicit Euler
end	
if timeMethod=="fracstep" then
	timeDisc = ThetaTimeStep(domainDisc,"FracStep")
end
if timeMethod=="alex" then
	timeDisc = ThetaTimeStep(domainDisc, "Alexander")
end

-- create the assembled operator for the solver
op = AssembledOperator(timeDisc)


op:init()

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

solverDesc =
{
	type = "newton",
	linSolver =
	{
		type = "bicgstab",
		precond =
		{
			type = "gmg",
			rap = false,
			smoother =
			{
				type = "ilu",
				beta = ilu_beta,
				consistentInterfaces = true
			},
			preSmooth = 1,
			postSmooth = 1,
			baseSolver = "lu",
			baseLevel = numPreRefs
		},
		convCheck =
		{
			type		= "standard",
			iterations	= 128,
			absolute	= 1e-12,
			reduction	= 0.5e-8,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		= 3,
		lambdaStart		= 1,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= 256,
		absolute	= 2.5e-12,
		reduction	= 1e-8,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Prepare the initial values
------------------------------------------------------------------------------------------

-- grid function for the solution
u = GridFunction (approxSpace)

-- set the initial guess
u:set (0)
Interpolate(0.0, u, "u")
Interpolate("ValueY", u, "v")
Interpolate("ValueP", u, "p")
Interpolate("ValueC", u, "c")

-- order the DoFs:
--OrderCuthillMcKee (approxSpace, true)
OrderLex (approxSpace,  "x")


------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0


vtk_file_name = file_name .. "-lev" .. numRefs
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end
	-- write start solution
	print("Writing start values")
	out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	out:select_nodal ("u,v", "velocity")
	out:select_nodal ("u", "u")
	out:select_nodal ("v", "v")
	out:select_nodal ("p", "p")
	out:select_nodal ("c", "c")


	out:print_subsets(vtk_file_name, u,allSubsets,0,0)
	
solver:init(op)


if solver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end
    
-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- choose time step
	do_dt = dt
		
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)
	
	-- prepare newton solver
	if solver:prepare(u) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 
	
	if solver:apply(u)  == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	
	end

	-- update new time
	time= timeDisc:future_time()
		
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecAssign(oldestSol, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	
	-- compute CFL number
	 
	CFL=cflNumber(u,do_dt)
	print("DT=" .. dt .. "");
	print("Time=" .. time .. "");
		
	
	if step % outputFactor == 0 then
	
		out = VTKOutput()
		out:clear_selection()
		out:select_all(false)
		out:select_nodal("u,v", "velocity")
		out:select_nodal("u", "u")
		out:select_nodal("v", "v")
		out:select_nodal("p", "p")
		out:select_nodal("c", "c")


		out:print_subsets(vtk_file_name, u,allSubsets,step,time)
		print(" ")
	end
	print("++++++ TIMESTEP " .. step .. "  END ++++++")
end

tAfter = os.clock()
solver:print_average_convergence()
print("Computation took " .. tAfter-tBefore .. " seconds.")

print("done.")
