-------------------------------------------------------------------------------------------------------
--
--  Stationary bubble
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "fvcr" or "stabil"
discType = util.GetParam("-type", "fvcr")
elemType = util.GetParam("-elem", "quads")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_symtri.ugx")
--		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_unstructured_tri_2x2.ugx")
	else		
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end
dt = util.GetParamNumber("-dt", 0.05)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 500)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",3)
La = util.GetParamNumber("-La",120)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. discType)
print("    dt           = " .. dt)
print("    numTimeSteps = " .. numTimeSteps)
print("    Laplace number = " .. La)

print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Problem Parameters
--------------------------------------------
--------------------------------------------

-- bubble radius
radius = 0.2

-- surface tension factor
sigma = 1

-- viscosity
mu = 0.1

-- density
rho = La * mu * mu / (2*radius * sigma)

-- gravitation
gravitation = 0

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
approxSpace:add_fct("p", "piecewise-constant") 

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

approxSpaceLevelSet = ApproximationSpace(dom)
approxSpaceLevelSet:add_fct("lsf", "Lagrange", 1)
approxSpaceLevelSet:add_fct("curvature", "piecewise-constant")

OrderLex(approxSpace, "lr");
OrderLex(approxSpaceLevelSet, "lr");

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", "fvcr")

-- set upwind
noUpwind = NavierStokesNoUpwind();
fullUpwind = NavierStokesFullUpwind();
weightedUpwind = NavierStokesWeightedUpwind(0.5);
NavierStokesDisc:set_conv_upwind(noUpwind)
	
NavierStokesDisc:set_peclet_blend(true)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_laplace(true)
NavierStokesDisc:set_stokes(true)

-- Level set class

lsDisc = FV1LevelSetDisc();

----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("Boundary")

----------------------------------
----------------------------------
-- Initial Level Set Function
----------------------------------
----------------------------------

function phiInit(x, y)
	return math.sqrt( (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5) )-0.2;
end

phiNew = GridFunction(approxSpaceLevelSet);
phiOld = GridFunction(approxSpaceLevelSet);
phiNew:set(0);
solfunctor = LuaUserNumber("phiInit");
Interpolate("phiInit", phiNew, "lsf");
VecAssign(phiOld,phiNew);

----------------------------------
----------------------------------
-- source 
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 0,0
end

rhs = LuaUserVector("source2d")

NavierStokesDisc:set_source(rhs)

density = LevelSetUserData(approxSpaceLevelSet,phiNew)
density:set_inside(rho)
density:set_outside(rho)
density:set_eval_type(0)

viscosity = LevelSetUserData(approxSpaceLevelSet,phiNew)
viscosity:set_inside(mu)
viscosity:set_outside(mu)

source = CRTwoPhaseSource(approxSpaceLevelSet,phiNew)
source:set_sigma(sigma)
source:set_density(density)
source:set_gravitation(gravitation)

NavierStokesDisc:set_source(source)
NavierStokesDisc:set_density(density)
NavierStokesDisc:set_kinematic_viscosity(viscosity)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(WallDisc)

-- create operator from discretization

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler
op = AssembledOperator(timeDisc)

op:init()

u = GridFunction(approxSpace)
u:set(0)

-- set level set velocity
velocity = GridFunctionVectorData(u, "u,v");
lsDisc:set_velocity(velocity)
lsDisc:set_dt(dt)
lsDisc:exact_curvature(1.0/radius)

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 0 end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

vanka = LineVanka(approxSpace)
vanka:set_num_steps(4,4,0,0,0,0)
-- vanka = Vanka()
vanka:set_damp(0.9)
-- vanka = Vanka()
-- vanka = BlockVanka(approxSpace)
-- vanka:set_blocksize(0.333334/2/2,0.333334/2/2,0.1)
-- vanka:update()
-- vanka:set_relax(0.9);

ILUT = CRILUT()
ILUT:set_threshold(1e-0,1e-1,1e-1,1e-1)
ILUT:set_info(true)
ILUT:set_damp(0.9)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(ILUT)
ilutSolver:set_convergence_check(ConvCheck(100000, 1e-12, 1e-14, false))

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(100000, 1e-12, 1e-14, false))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-12)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(true)

vankaBase = LinearSolver()
vankaBase:set_preconditioner(Vanka())
vankaBase:set_convergence_check(baseConvCheck)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(vankaSolver)
gmg:set_smoother(ILUT)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(0.8)
-- gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(vanka)
-- BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-12, 1e-14, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-12, 1e-14, true))

-- choose a solver
solver = BiCGStabSolver
solver = vankaSolver
solver = gmgSolver
solver = ilutSolver

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10000)
newtonConvCheck:set_minimum_defect(1e-12)
newtonConvCheck:set_reduction(1e-14)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.7)
newtonLineSearch:set_accept_best(false)

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

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("twoPhase", u,0,0)

lsfilename = "lsf"
outls = VTKOutput()
step=0;
time=0;
outls:select_element("lsf", "lsf")
outls:print(lsfilename, phiNew, step, time)

-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

	-- store grid function in vector of  old solutions
	solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(uOld, time)
	
	lsDisc:compute_curvature(phiNew);

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
		
		lsDisc:advect_lsf(phiNew,phiOld);
		lsDisc:compute_curvature(phiNew);
		VecAssign(phiOld,phiNew);
		
		-- plot solution

		out:clear_selection()
		out:print("twoPhase", u,step,time)
		outls:print(lsfilename, phiNew, step, time)
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
SaveVectorForConnectionViewer(u, "u.vec")

print("done.")
