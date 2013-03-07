-------------------------------------------------------------------------------------------------------
--
--  Cylinder Problem using staggered grid discretization on Crouzeix-Raviart type elements
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then gridName = util.GetParam("-grid", "grids/cylinder.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

requiredSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall", "CylinderWall"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
approxSpace:add_fct("p", "piecewise-constant") 

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

u = GridFunction(approxSpace)

tBefore = os.clock()
-- OrderCRCuthillMcKee(approxSpace,u,true)
-- CROrderCuthillMcKee(approxSpace,u,true,false,false,true)
-- CROrderSloan(approxSpace,u,false,false,true)
-- CROrderKing(approxSpace,u,true,false,false,true)
CROrderMinimumDegree(approxSpace,u,true,false,true)
-- OrderLex(approxSpace, "lr");
tAfter = os.clock()
print("Ordering took " .. tAfter-tBefore .. " seconds.");

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", "staggered")
	
-- set upwind
noUpwind = NavierStokesNoUpwind();
fullUpwind = NavierStokesFullUpwind();
weightedUpwind = NavierStokesWeightedUpwind(0.6);
NavierStokesDisc:set_conv_upwind(weightedUpwind)

NavierStokesDisc:set_peclet_blend(false)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_stokes(false)
NavierStokesDisc:set_laplace(true)
NavierStokesDisc:set_kinematic_viscosity(1.0e-3);

----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function inletVel2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return 4 * Um * y * (H-y) / (H*H), 0.0
end

OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDisc:add("Outlet")

-- setup inlet
InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Inlet")

WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("UpperWall,LowerWall,CylinderWall")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 0,0
end

rhs = LuaUserVector("source2d")

NavierStokesDisc:set_source(rhs)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)
domainDisc:add(WallDisc)
domainDisc:add(OutletDisc)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------

-- create operator from discretization
op = AssembledOperator(domainDisc)
op:init()

u:set(0)

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 0 end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

-- OrderLex(approxSpace, "lr");

vanka = LineVanka(approxSpace)
vanka:set_num_steps(2,2,2,2)
vanka:set_damp(0.9)

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-1, true))

CRILUT = CRILUT()
CRILUT:set_threshold(1e-0,1e-1,1e-1,1e-1)
-- CRILUT:set_threshold(1e-0)
CRILUT:set_damp(0.9)
CRILUT:set_info(true)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(CRILUT)
ilutSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-1, true))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-7)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(false)

baseVanka = Vanka()
baseVanka:set_relax(0.9)
vankaBase = LinearSolver()
vankaBase:set_preconditioner(baseVanka)
vankaBase:set_convergence_check(baseConvCheck)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(vankaBase)
gmg:set_smoother(CRILUT)
gmg:set_cycle_type(1)
gmg:set_damp(MinimalResiduumDamping())
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
-- gmg:add_restriction_post_process(AverageComponent(approxSpace,"p"))
-- gmg:add_prolongation_post_process(AverageComponent(approxSpace,"p"))

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-1, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-1, true))

-- choose a solver
solver = gmgSolver
-- solver = ilutSolver
-- solver = vankaSolver

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(50)
newtonConvCheck:set_minimum_defect(1e-5)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(false)

dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed."); exit();
end

SaveVectorForConnectionViewer(u, "StartSolution.vec")

tBefore = os.clock()

if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("CRCylinderSolution", u)

print("done.")

