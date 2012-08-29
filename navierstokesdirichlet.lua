-------------------------------------------------------------------------------------------------------
--
--  Lua - Script to test the Navier-Stokes implementation
--
--  Author: Christian Wehner
--
--	A theoretical example to test the Navier-Stokes discretization.
--	The boundary conditions are inflow boundary conditions (Dirichlet conditions for the velocity)
--  on the whole boundary.
--  The analytical solution can be constructed e.g. via maple:
--
--  # construct u and v solution so that velocity is divergence-free
--  g:=x*x*x+x*y+y*y;
--  u:=diff(g,y);
--  v:=-diff(g,x); 
--  # chose p
--  p:=x;
--  # rhs is chosen so that Navier-Stokes system is fulfilled
--  rhsu:=factor(simplify(-diff(diff(u,x),x)-diff(diff(u,y),y))+u*diff(u,x)+v*diff(u,y)+diff(p,x));
--  rhsv:=factor(simplify(-diff(diff(v,x),x)-diff(diff(v,y),y))+u*diff(v,x)+v*diff(v,y)+diff(p,y));
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "staggered" or "stabil"
discType = util.GetParam("-type", "staggered")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
-- gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_1x1.ugx")
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

elemDisc = NavierStokes(fctUsed, "Inner")

if discType=="staggered" then
	
	elemDisc:set_disc_scheme("staggered");
	-- set upwind
	noUpwind = NavierStokesCRNoUpwind();
	fullUpwind = NavierStokesCRFullUpwind();
	elemDisc:set_conv_upwind(noUpwind)
	
else

    -- stabilization scheme
	elemDisc:set_disc_scheme("stab");
	
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
	elemDisc:set_stabilization(stab)
	
	-- set upwind
	elemDisc:set_conv_upwind(upwind)

end

elemDisc:set_peclet_blend(false)
elemDisc:set_exact_jacobian(false)
elemDisc:set_stokes(false)
elemDisc:set_laplace(true)
elemDisc:set_kinematic_viscosity(1.0);


----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
	return x+2*y		
end

function vsol2d(x,y,t)
	return -4*x*x*x-y
end

function psol2d(x,y,t)
	return x*y
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

uSolution = LuaUserNumber("usol"..dim.."d")
vSolution = LuaUserNumber("vsol"..dim.."d")
pSolution = LuaUserNumber("psol"..dim.."d")

if discType=="stabil" then
	LuaInletDisc = NavierStokesInflow("u,v,p", "Inner")
else
	LuaInletDisc = CRNavierStokesInflow("u,v", "Inner")
end
LuaInletDisc:add("inletVel"..dim.."d", "Boundary")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return x+y-8*x*x*x,25*x-8*x*x*x-24*x*x*y+y
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
vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(10000, 1e-10, 1e-12, true))

-- create Linear Solver
linSolver = BiCGStab()
if type=="stabil" then
	linSolver:set_preconditioner(gmg)
else
	linSolver:set_preconditioner(vanka)
end
linSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-12, true))

-- choose a solver
-- solver = linSolver
solver = vankaSolver
-- solver = LU()
-- solver = linSolver

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(20)
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
out:print("NavierStokesSolution", u)

print("done.")
