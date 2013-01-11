------------------------------------------------------------------
--
--   Lua - Script to test the Navier-Stokes implementation
--
--	This script sets up a problem for the Navier-Stokes discretization
--	on staggered grid using Crouzeix-Raviart type elements
--
--   Author: Christian Wehner
--
-------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "staggered" or "stabil"
discType = util.GetParam("-type", "staggered")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_1x1.ugx")
-- gridName = util.GetParam("-grid", "unit_square/unit_square_unstructured_tris_coarse.ugx")
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. discType)
print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

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

NavierStokesDisc = NavierStokes(fctUsed, "Inner", discType)

if discType=="staggered" then
	-- set upwind
	noUpwind = NavierStokesNoUpwind();
	fullUpwind = NavierStokesFullUpwind();
	NavierStokesDisc:set_conv_upwind(noUpwind)
	
else
	--upwind = NavierStokesNoUpwind();
	--upwind = NavierStokesFullUpwind();
	--upwind = NavierStokesSkewedUpwind();
	upwind = NavierStokesLinearProfileSkewedUpwind();
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
	NavierStokesDisc:set_stabilization(stab)
	
	-- set upwind
	NavierStokesDisc:set_conv_upwind(upwind)

end

NavierStokesDisc:set_peclet_blend(false)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_stokes(true)
NavierStokesDisc:set_laplace(true)
NavierStokesDisc:set_kinematic_viscosity(1.0);


-----------------------------------
----------------------------------
-- Dirichlet boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
	  return 0
--    return x
--    return 2*x*y*(x-1)*(y-1)*(-x*(x-1)*(2*y-1))
--    return y*y
--    return y
end

function vsol2d(x,y,t)
--    return -3*x*x
--    return -3*x*x-y
--		return 2*x*y*(x-1)*(y-1)*(y*(y-1)*(2*x-1))
		return x*x
--    return x
end

function psol2d(x,y,t)
    return 0
end

uSolution = LuaUserNumber("usol"..dim.."d")
vSolution = LuaUserNumber("vsol"..dim.."d")
pSolution = LuaUserNumber("psol"..dim.."d")

dirichletBnd = DirichletBoundary()
dirichletBnd:add(uSolution, "u", "Boundary")
dirichletBnd:add(vSolution, "v", "Boundary")
dirichletBnd:add(pSolution, "p", "Boundary")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	  return 0,-2
--    return 4*(2*y-1)*(-6*x*y*y+6*x*x*y*y+y*y+6*x*y-6*x*x*y-y+3*x*x*x*x+3*x*x-6*x*x*x),-4*(2*x-1)*(6*x*x*y*y-6*x*x*y+x*x-6*x*y*y+6*x*y-x-6*y*y*y+3*y*y+3*y*y*y*y)
--    return -2,-2
--    return  0,0
end

rhs = LuaUserVector("source2d")

NavierStokesDisc:set_source(rhs)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(dirichletBnd)

op = AssembledOperator(domainDisc)

u = GridFunction(approxSpace)
u:set(0)

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 0 end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(LU())
gmg:set_smoother(ILU())
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
--gmg:set_debug(dbgWriter)

vanka = Vanka()
-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(10000, 1e-10, 1e-12, true))

-- create Linear Solver
linSolver = BiCGStab()
if discType=="stabil" then
	linSolver:set_preconditioner(gmg)
else
	linSolver:set_preconditioner(vanka)
end
linSolver:set_convergence_check(ConvCheck(100000, 1e-10, 1e-12, true))

-- choose a solver
if discType=="stabil" then
	solver = linSolver
else
	solver = vankaSolver
end

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(5e-10)
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
newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed."); exit();
end

SaveVectorForConnectionViewer(u, "StartSolution.vec")

if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

l2error = L2Error(uSolution, u, "u", 0.0, 1, "Inner")
write("L2Error in u component is "..l2error .."\n");
l2error = L2Error(vSolution, u, "v", 0.0, 1, "Inner")
write("L2Error in v component is "..l2error .."\n");
l2error = L2Error(pSolution, u, "p", 0.0, 1, "Inner")
write("L2Error in p component is "..l2error .."\n");


out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("StokesSolution", u)

print("done.")
