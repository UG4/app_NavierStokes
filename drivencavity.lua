-------------------------------------------------------------------------------------------------------
--
--  Driven cavity problem
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "staggered" or "stabil"
discType = util.GetParam("-type", "staggered")
boolStat = util.GetParamNumber("-stat", 1)
elemType = util.GetParam("-elem", "quad")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "grids/dc_tri.ugx")
	else
		gridName = util.GetParam("-grid", "grids/dc_quads.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end
dt = util.GetParamNumber("-dt", 0.1)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",3)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. discType)
if boolStat==false then
	print("    dt           = " .. dt)
	print("    numTimeSteps = " .. numTimeSteps)
end
print("    grid       	= " .. gridName)
print("    stationary  	= " .. boolStat)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

if discType=="staggered" then
	if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
	if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
	if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
	approxSpace:add_fct("p", "piecewise-constant") 
else
	if dim >= 1 then approxSpace:add_fct("u", "Lagrange", 1) end
	if dim >= 2 then approxSpace:add_fct("v", "Lagrange", 1) end
	if dim >= 3 then approxSpace:add_fct("w", "Lagrange", 1) end
	approxSpace:add_fct("p", "Lagrange", 1)
end

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

elemDisc = NavierStokes(fctUsed, "Inner")

if discType=="staggered" then
	
	elemDisc:set_disc_scheme("staggered");
	-- set upwind
	noUpwind = NavierStokesCRNoUpwind();
	fullUpwind = NavierStokesCRFullUpwind();
	weightedUpwind = NavierStokesCRWeightedUpwind(0.5);
	elemDisc:set_conv_upwind(weightedUpwind)
	
else

    -- stabilization scheme
	elemDisc:set_disc_scheme("stab");
	
	--upwind = NavierStokesNoUpwind();
	upwind = NavierStokesFullUpwind();
	--upwind = NavierStokesSkewedUpwind();
	--upwind = NavierStokesLinearProfileSkewedUpwind();
	--upwind = NavierStokesRegularUpwind();
	--upwind = NavierStokesPositiveUpwind();
	-- chose stabilization 
	stab = NavierStokesFIELDSStabilization()
	--stab = NavierStokesFLOWStabilization()
	-- ... and set the upwind
	stab:set_upwind(upwind)
	--stab:set_diffusion_length("RAW")
	stab:set_diffusion_length("FIVEPOINT")
	--stab:set_diffusion_length("COR")
	
	-- set stabilization
	elemDisc:set_stabilization(stab)
	
	-- set upwind
	elemDisc:set_conv_upwind(upwind)

end

elemDisc:set_peclet_blend(true)
elemDisc:set_exact_jacobian(false)
elemDisc:set_stokes(false)
elemDisc:set_laplace(true)
elemDisc:set_kinematic_viscosity(1.0/3200.0);


----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
	return 1		
end

function vsol2d(x,y,t)
	return 0
end

function psol2d(x,y,t)
	return 0
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

uSolution = LuaUserNumber("usol"..dim.."d")
vSolution = LuaUserNumber("vsol"..dim.."d")
pSolution = LuaUserNumber("psol"..dim.."d")

if discType=="stabil" then
	LuaInletDisc = NavierStokesInflow("u,v,p", "Inner")
	WallDisc = NavierStokesWall("u,v,p")
else
	LuaInletDisc = CRNavierStokesInflow("u,v", "Inner")
	WallDisc = CRNavierStokesWall("u,v,p")
end
LuaInletDisc:add("inletVel"..dim.."d", "Top")
WallDisc:add("Left,Right,Bottom")


----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 0,0
end

rhs = LuaUserVector("source2d")

elemDisc:set_source(rhs)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(LuaInletDisc)
domainDisc:add(WallDisc)

-- create operator from discretization

if boolStat==1 then
	-- operator is stationary
	op = AssembledOperator(domainDisc)
else
	-- create time discretization
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1.0) -- 1.0 is implicit euler
	op = AssembledOperator(timeDisc)
end
op:init()

u = GridFunction(approxSpace)
u:set(0)

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 0 end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

vanka = Vanka()
vanka:set_damp(0.95)

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(StandardConvergenceCheck(100000, 1e-7, 2.5e-1, true))

baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-7)
baseConvCheck:set_reduction(1e-2)
baseConvCheck:set_verbose(false)

vankaBase = LinearSolver()
vankaBase:set_preconditioner(Vanka())
vankaBase:set_convergence_check(baseConvCheck)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(vankaBase)
gmg:set_smoother(vanka)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(8*numRefs)
gmg:set_num_postsmooth(8*numRefs)
-- gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(0.8)
-- gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(StandardConvergenceCheck(100000, 1e-7, 2.5e-1, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(StandardConvergenceCheck(10000, 1e-7, 2.5e-1, true))

-- choose a solver
solver = BiCGStabSolver
solver = gmgSolver

newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10000)
newtonConvCheck:set_minimum_defect(1e-6)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)

dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
-- newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed."); exit();
end

SaveVectorForConnectionViewer(u, "StartSolution.vec")

-- if newtonSolver:apply(u) == false then
--	 print ("Newton solver apply failed."); exit();
-- end

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- start
time = 0.0
step = 0

-- setup the lua functions ...
function Pressure_StartValue2d(x, y, t) return 0.0 end
function VelX_StartValue2d(x, y, t) return 0.0 end
function VelY_StartValue2d(x, y, t)	return 0.0 end

-- Now interpolate the function
time = 0.0
Interpolate("Pressure_StartValue"..dim.."d", u, "p", time);
Interpolate("VelX_StartValue"..dim.."d", u, "u", time);
Interpolate("VelY_StartValue"..dim.."d", u, "v", time);

-- filename
filename = "Sol"

-- write start solution
-- print("Writing start values")
-- out = VTKOutput()
-- out:print(filename, u, step, time)

-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()
if boolStat==1 then

	newtonSolver:init(op)
		if newtonSolver:prepare(u) == false then
		print ("Newton solver prepare failed."); exit();
	end

	SaveVectorForConnectionViewer(u, "StartSolution.vec")

	if newtonSolver:apply(u) == false then
		print ("Newton solver apply failed."); exit();
	end

else

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
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
	
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 

		-- update new time
		time = solTimeSeries:time(0) + do_dt
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()

		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
	
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)

		print("++++++ TIMESTEP " .. step .. "  END ++++++");
		
		-- plot solution

		out = VTKOutput()
		out:clear_selection()
		out:select_all(false)
		out:select_element("u,v", "velocity")
		out:select_element("u", "u")
		out:select_element("v", "v")
		out:select_element("p", "p")
		out:print("timeDCSolution", u)

	end

end

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

-- plot solution

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("DCSolution", u)

print("done.")