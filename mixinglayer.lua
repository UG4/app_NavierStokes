-------------------------------------------------------------------------------------------------------
--
--  Mixing layer problem (test case for turbulence)
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "staggered" or "stabil"
discType = util.GetParam("-type", "staggered")
elemType = util.GetParam("-elem", "quads")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then 
--		gridName = util.GetParam("-grid", "grids/dc_tri.ugx")
	else
		gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end
dt = util.GetParamNumber("-dt", 0.1)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 100)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",3)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    dt           = " .. dt)
print("    numTimeSteps = " .. numTimeSteps)
print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
IdentifySubsets(dom,"Left","Right");

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

approxSpaceVorticity = ApproximationSpace(dom)
approxSpaceVorticity:add_fct("c", "Crouzeix-Raviart", 1)

vort = GridFunction(approxSpaceVorticity);

u = GridFunction(approxSpace);
-- OrderCRCuthillMcKee(approxSpace,u,true);
-- OrderLex(approxSpace, "lr");


-- Setup from John LES book computation on [-1 1] square
sigma0=1.0/14.0
winf=1
cnoise=0.01
viscosity=140000

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

elemDisc = NavierStokes(fctUsed, "Inner", "staggered")

-- set upwind
noUpwind = NavierStokesCRNoUpwind();
fullUpwind = NavierStokesCRFullUpwind();
weightedUpwind = NavierStokesCRWeightedUpwind(0.5);
elemDisc:set_conv_upwind(fullUpwind)	

elemDisc:set_peclet_blend(true)
elemDisc:set_exact_jacobian(false)
elemDisc:set_stokes(false)
elemDisc:set_laplace(false)

----------------------------------
----------------------------------
-- Viscosity Data
----------------------------------
----------------------------------

viscosityData = CRSmagorinskyTurbViscData(approxSpace,u,0.1)
viscosityData:set_kinematic_viscosity(1/140000);
elemDisc:set_kinematic_viscosity(viscosityData);
elemDisc:set_laplace(false);


----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

-- OutletDiscTop = CRNavierStokesNoNormalStressOutflow(elemDisc)
OutletDiscTop = CRNavierStokesSymBC(elemDisc)
OutletDiscTop:add("Top")
-- OutletDiscBottom = CRNavierStokesNoNormalStressOutflow(elemDisc)
OutletDiscBottom = CRNavierStokesSymBC(elemDisc)
OutletDiscBottom:add("Bottom")

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
domainDisc:add(OutletDiscTop)
domainDisc:add(OutletDiscBottom)

-- create operator from discretization

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(0.5) -- 1.0 is implicit euler
op = AssembledOperator(timeDisc)

op:init()

u:set(0)

function StartValue_u2d(x,y,t) 
	return winf*math.tanh(2*y/sigma0)+cnoise*winf*(-8*y*math.exp(-(2*y/sigma0)*(2*y/sigma0))*(math.cos(8*math.pi*x)+math.cos(20*math.pi*x)))/sigma0/sigma0;
end

function StartValue_v2d(x,y,t) 
	return cnoise*winf*(-math.exp(-4*y*y/sigma0/sigma0)*(-8*math.sin(8*math.pi*x)*math.pi-20*math.sin(20*math.pi*x)*math.pi))
end

function StartValue_p2d(x,y,t) return 0 end

Interpolate("StartValue_u2d", u, "u")
Interpolate("StartValue_v2d", u, "v")
Interpolate("StartValue_p2d", u, "p")

-- vanka = LineVanka(approxSpace)
-- vanka:set_num_steps(4,4,0,0,0,0)
-- vanka = CRILU()
--vanka = Vanka()
-- vanka:set_damp(0.9)
vanka = Vanka()
-- vanka = BlockVanka(approxSpace)
-- vanka:set_blocksize(0.333334/2/2,0.333334/2/2,0.1)
-- vanka:update()
-- vanka:set_relax(0.9);
vanka:set_damp(0.9)

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, true))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-5)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(false)

CRILUT = CRILUT()
CRILUT:set_threshold(1e-0,1e-2,1e-2,1e-2)
CRILUT:set_damp(1)
-- CRILUT:set_info(true)
ILUTBase = ILUT()
ILUTBase:set_threshold(0)
ILUTBase:set_info(false)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(CRILUT)
ilutSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, true))

vankaBase = LinearSolver()
vankaBase:set_preconditioner(Vanka())
vankaBase:set_convergence_check(baseConvCheck)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(ilutSolver)
gmg:set_smoother(vanka)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
-- gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(0.8)
-- gmg:set_damp(MinimalEnergyDamping())

-- gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(vanka)
-- BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-5, 1e-1, true))

-- choose a solver
solver = BiCGStabSolver
solver = vankaSolver
-- solver = gmgSolver
solver = ilutSolver

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10000)
newtonConvCheck:set_minimum_defect(1e-5)
newtonConvCheck:set_reduction(1e-10)
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
-- newtonSolver:set_line_search(newtonLineSearch)
-- newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)
newtonSolver:add_step_update(viscosityData)

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

-- Now interpolate the function
time = 0.0

-- filename
filename = "Sol"

-- compute initial vorticity
vorticity(vort,u)

-- write start solution
print("Writing start values")
out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("MixingLayer", u,0,0)

outv = VTKOutput()
outv:select_element("c","c")
outv:print("vorticity", vort,0,0)

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
	
	-- compute vorticity
	vort:set(0)
	vorticity(vort,u)
	
	-- plot solution

	out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	out:select_element("u,v", "velocity")
	out:select_element("u", "u")
	out:select_element("v", "v")
	out:select_element("p", "p")
	out:print("MixingLayer", u,step,time)
	outv = VTKOutput()
	outv:select_element("c","c")
	outv:print("vorticity", vort,step,time)

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
