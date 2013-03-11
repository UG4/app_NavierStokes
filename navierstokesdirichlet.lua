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
--  rhsu:=factor(simplify(-1/R*diff(diff(u,x),x)-1/R*diff(diff(u,y),y)+u*diff(u,x)+v*diff(u,y)+diff(p,x)));
--  rhsv:=factor(simplify(-1/R*diff(diff(v,x),x)-1/R*diff(diff(v,y),y)+u*diff(v,x)+v*diff(v,y)+diff(p,y)));
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "fvcr" or "stabil"
discType = util.GetParam("-type", "fvcr")
order = util.GetParamNumber("-order", 1)
vorder = util.GetParamNumber("-vorder", order)
porder = util.GetParamNumber("-porder", order-1)

InitUG(dim, AlgebraType("CPU", 1));

R=100

if 	dim == 2 then
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_1x1.ugx")
-- gridName = util.GetParam("-grid", "grids/stretched_quads.ugx") 
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. discType)
print("    grid       	= " .. gridName)
print("    v ansatz order = " ..vorder)
print("    p ansatz order = " ..porder)

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

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"}; FctCmp = {"u", "p"};
elseif  dim == 2 then VelCmp = {"u", "v"}; FctCmp = {"u", "v", "p"};
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"};
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

if discType=="fvcr" then
	if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
	if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
	if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
	approxSpace:add_fct("p", "piecewise-constant") 
end
if discType=="fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
end
if discType=="fv" or discType=="fe" then
	approxSpace:add_fct(VelCmp, "Lagrange", vorder)
	if porder==0 then
		approxSpace:add_fct("p", "piecewise-constant") 
    else 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
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

if discType=="fvcr" then
	-- set upwind
	noUpwind = NavierStokesNoUpwind();
	fullUpwind = NavierStokesFullUpwind();
	NavierStokesDisc:set_conv_upwind(fullUpwind)
	NavierStokesDisc:set_defect_upwind(false);
	NavierStokesDisc:set_defect_upwind(true);
end	
if discType=="fv1" then
	upwind = NavierStokesNoUpwind();
	-- upwind = NavierStokesFullUpwind();
	-- upwind = NavierStokesSkewedUpwind();
	-- upwind = NavierStokesLinearProfileSkewedUpwind();
	--upwind = NavierStokesRegularUpwind();
	--upwind = NavierStokesPositiveUpwind();
	-- chose stabilization 
	stab = NavierStokesFIELDSStabilization()
	stab = NavierStokesFLOWStabilization()
	-- ... and set the upwind
	stab:set_upwind(upwind)
	--stab:set_diffusion_length("RAW")
	--stab:set_diffusion_length("FIVEPOINT")
	stab:set_diffusion_length("COR")
	
	-- set stabilization
	NavierStokesDisc:set_stabilization(stab)
	
	-- set upwind
	NavierStokesDisc:set_conv_upwind(upwind)

end

NavierStokesDisc:set_peclet_blend(false)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_stokes(false)
NavierStokesDisc:set_laplace(true)
NavierStokesDisc:set_kinematic_viscosity(1.0/R);



----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
	return 2*x^2*(1-x)^2*(y*(1-y)^2-y^2*(1-y))		
--	return y;
end

function vsol2d(x,y,t)
	return -2*y^2*(1-y)^2*(x*(1-x)^2-x^2*(1-x))
--	return -x
end

function psol2d(x,y,t)
--	return x^3+y^3+0.5
	return 0
--	return x+y
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

uSolution = LuaUserNumber("usol"..dim.."d")
vSolution = LuaUserNumber("vsol"..dim.."d")
pSolution = LuaUserNumber("psol"..dim.."d")

if discType=="stabil" then
	InletDisc = NavierStokesInflow(NavierStokesDisc)
else
	InletDisc = NavierStokesInflow(NavierStokesDisc)
end
InletDisc:add("inletVel"..dim.."d", "Boundary")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 
--	3*x*x+ -- dx p term
(-4*y+12*x^2+12*y^2-24*x^3-8*y^3+24*x*y-72*x*y^2+48*x*y^3-48*x^2*y+72*x^2*y^2-48*x^2*y^3+48*x^3*y-24*x^4*y+12*x^4-16*x^3*y^3*R-20*x^4*y^2*R+36*x^5*y^2*R+28*x^3*y^4*R-28*x^6*y^2*R-24*x^3*y^5*R+80*x^4*y^3*R-140*x^4*y^4*R+120*x^4*y^5*R-144*x^5*y^3*R+252*x^5*y^4*R-216*x^5*y^5*R+112*x^6*y^3*R-32*x^7*y^3*R+8*x^7*y^2*R-196*x^6*y^4*R+168*x^6*y^5*R-40*x^4*y^6*R+72*x^5*y^6*R-56*x^6*y^6*R+56*x^7*y^4*R-48*x^7*y^5*R+16*x^7*y^6*R+8*x^3*y^6*R+4*x^3*y^2*R)/R
	,
--	3*y*y+ -- dy p term
(4*x-12*x^2-12*y^2+8*x^3+24*y^3-24*x*y+48*x*y^2-48*x*y^3+72*x^2*y-72*x^2*y^2-48*x^3*y+48*x^3*y^2+24*x*y^4-12*y^4-16*x^3*y^3*R+80*x^3*y^4*R-144*x^3*y^5*R+28*x^4*y^3*R-140*x^4*y^4*R+252*x^4*y^5*R-24*x^5*y^3*R+120*x^5*y^4*R-216*x^5*y^5*R+8*x^6*y^3*R-40*x^6*y^4*R+72*x^6*y^5*R-196*x^4*y^6*R+168*x^5*y^6*R-56*x^6*y^6*R+112*x^3*y^6*R-20*x^2*y^4*R+36*x^2*y^5*R-28*x^2*y^6*R-32*x^3*y^7*R+56*x^4*y^7*R-48*x^5*y^7*R+16*x^6*y^7*R+8*x^2*y^7*R+4*x^2*y^3*R)/R
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
domainDisc:add(InletDisc)

op = AssembledOperator(domainDisc)

u = GridFunction(approxSpace)
u:set(0)

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 1 end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

vanka = Vanka()
vanka:set_damp(1)

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(10000, 1e-9, 1e-12, true))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-7)
baseConvCheck:set_reduction(1e-3)
baseConvCheck:set_verbose(false)

vankaBase = LinearSolver()
vankaBase:set_preconditioner(DiagVanka())
vankaBase:set_convergence_check(baseConvCheck)

crilut = CRILUT()
crilut:set_threshold(1e-1,1e-4,1e-4,1e-4)
crilut:set_damp(1)
-- crilut:set_info(true)
ilutBase = ILUT()
ilutBase:set_threshold(1e-4)
ilutBase:set_info(false)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(ilutBase)
ilutSolver:set_convergence_check(ConvCheck(100000, 5e-8, 1e-1, false))


gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(ilutSolver)
if discType=="fv1" then
	ilut = ILUT()
	ilut:set_threshold(1e-3)
	ilut:set_info(false)
	gmg:set_smoother(ilut)
end
if discType=="fvcr" then
	gmg:set_smoother(vanka)
end
gmg:set_cycle_type(1)
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:add_prolongation_post_process(AverageComponent(approxSpace,"p"))
-- gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-3, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-2, true))


-- choose a solver
if discType=="fvcr" or discType=="fv1" then
	solver = gmgSolver
end
if discType=="fe" or discType=="fv" then
	if porder==0 then
		solver = vankaSolver
	else
		solver = LU()
	end
end

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(100)
newtonConvCheck:set_minimum_defect(1e-7)
if discType=="fv" or discType=="fe" then
	newtonConvCheck:set_reduction(0.999)
else 
	newtonConvCheck:set_reduction(1e-10)
end
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(30)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.85)
newtonLineSearch:set_accept_best(true)

dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

tBefore = os.clock()

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

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");
