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
-- 
-- u:=sin(2*pi*(x+t))*cos(2*pi*y);
-- v:=-cos(2*pi*(x+t))*sin(2*pi*y);
-- # chose p
-- p:=0.25*x*x;
-- # rhs is chosen so that Navier-Stokes system is fulfilled
-- rhsu:=factor(diff(u,t)+simplify(-1/R*diff(diff(u,x),x)-1/R*diff(diff(u,y),y)+u*diff(u,x)+v*diff(u,y)+diff(p,x)));
-- rhsv:=factor(diff(v,t)+simplify(-1/R*diff(diff(v,x),x)-1/R*diff(diff(v,y),y)+u*diff(v,x)+v*diff(v,y)+diff(p,y)));
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "staggered" or "stabil"
discType = util.GetParam("-type", "staggered")

InitUG(dim, AlgebraType("CPU", 1));

dim = util.GetParamNumber("-dim", 2)
type = util.GetParam("-type", "fvcr")
order = util.GetParamNumber("-order", 1)
vorder = util.GetParamNumber("-vorder", order)
porder = util.GetParamNumber("-porder", vorder-1)
bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "no", "Upwind type")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
adaptive    = util.HasParamOption("-adaptive", "If defined, adaptive grid is used")
numAdRefs   = util.GetParamNumber("-numAdaptRefs", 2)
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "COR", "Diffusion length type")
linred      = util.GetParam("-linred", 1e-1 , "Linear reduction")
nlintol     = util.GetParam("-nlintol", 1e-5, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")
elemType = util.GetParam("-elem", "quad")
dt = util.GetParamNumber("-dt", 0.05)
R = util.GetParamNumber("-R",100)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 50)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",4)

print(" Chosen Parameters:")
print("    dim                 = " .. dim)
print("    numTotalRefs        = " .. numRefs)
print("    numPreRefs          = " .. numPreRefs)
if adaptive==true then
	print("    numAdRefs           = " .. numAdRefs)
end
print("    dt                  = " .. dt)
print("    numTimeSteps        = " .. numTimeSteps)
print("    type                = " .. type)
print("    adaptive            = " .. tostring(adaptive))
print("    v ansatz order      = " ..vorder)
print("    p ansatz order      = " ..porder)
print("    no laplace          = " .. tostring(bNoLaplace))
print("    exact jacobian      = " .. tostring(bExactJac))
print("    peclet blend        = " .. tostring(bPecletBlend))
print("    upwind              = " .. upwind)
if type=="fv1" then
	print("    pac upwind          = " .. tostring(bPac))
	print("    stab                = " .. stab)
	print("    diffLength          = " .. diffLength)
end
print("    linear reduction    = " .. linred)
print("    linear tolerance    = " .. lintol)
print("    nonlinear reduction = " .. nlinred)
print("    nonlinear tolerance = " .. nlintol)
print("    Reynolds-nr         = " .. R)


if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "unit_square/unit_square_tri_4bnd.ugx")
	else
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2_four_bnd.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end


--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

approxSpace = ApproximationSpace(dom)

if 		dim == 1 then VelCmp = {"u"} FctCmp = {"u", "p"}
elseif  dim == 2 then VelCmp = {"u", "v"} FctCmp = {"u", "v", "p"}
elseif  dim == 3 then VelCmp = {"u", "v", "w"} FctCmp = {"u", "v", "w", "p"}
else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end

if type == "fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
elseif type == "fv" then
	approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
	approxSpace:add_fct("p", "Lagrange", porder) 
elseif type == "fe" then
	if porder==0 then
		if vorder==1 then
			approxSpace:add_fct(VelCmp, "Crouzeix-Raviart",1)
		else
			approxSpace:add_fct(VelCmp, "Lagrange", vorder)
		end
		approxSpace:add_fct("p", "piecewise-constant") 
	else
		approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
elseif type=="fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else print("Disc Type '"..type.."' not supported."); exit(); end

approxSpace:init_levels()

-- finally we print some statistic on the distributed dofs
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, type)

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(1.0/R)




--upwind if available
if type == "fv1" or type == "fvcr" then
	NavierStokesDisc:set_upwind(upwind)
	NavierStokesDisc:set_peclet_blend(bPecletBlend)
end

-- fv1 must be stablilized
if type == "fv1" then
	NavierStokesDisc:set_stabilization(stab, diffLength)
--  use PAC upwind or not
	NavierStokesDisc:set_pac_upwind(bPac)
end

-- fe must be stabilized for (Pk, Pk) space
if type == "fe" and porder == vorder then
	NavierStokesDisc:set_stabilization(3)
end

----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
	return math.sin(2*math.pi*(x+t))*math.cos(2*math.pi*y)		
end

function vsol2d(x,y,t)
	return -math.cos(2*math.pi*(x+t))*math.sin(2*math.pi*y)
end

function psol2d(x,y,t)
	return math.sin(x)*math.sin(y)*math.sin(t)
end

function inletVel2d(x, y, t)
	return -y*math.cos(t),x*math.cos(t)
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Top,Bottom,Left,Right")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 
	2*math.pi*(math.cos(2*math.pi*(x+t))*math.cos(2*math.pi*y)*R+4*math.sin(2*math.pi*(x+t))*math.pi*math.cos(2*math.pi*y)+math.sin(2*math.pi*(x+t))*math.cos(2*math.pi*(x+t))*R)/R
	,
	-2*math.pi*math.sin(2*math.pi*y)*(-math.sin(2*math.pi*(x+t))*R+4*math.cos(2*math.pi*(x+t))*math.pi-math.cos(2*math.pi*y)*R)/R
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

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(0.5) -- 1.0 is implicit euler
op = AssembledOperator(timeDisc)

op:init()

u = GridFunction(approxSpace)
if adaptive==true then
	for i=0, dim do
		u:add_transfer(P1LocalTransfer(i))
	end
end

function StartValue_u(x,y,t) return 0 end
function StartValue_v(x,y,t) return 0 end
function StartValue_p(x,y,t) return 0 end

Interpolate("usol2d", u, "u")
Interpolate("vsol2d", u, "v")
Interpolate("psol2d", u, "p")

-- first refinement

if adaptive==true then
	refiner = AdaptiveRegularDomainRefiner(dom)
	for i=1,numAdRefs do
		MarkForAdaption_GradientIndicator(refiner, u, "u", 1e-5, 0.75, 0.75, 100)
		refiner:refine() 
		refiner:clear_marks()
	end
end

-- finally we print some statistic on the distributed dofs
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

-- Interpolate("usol2d", u, "u")
-- Interpolate("vsol2d", u, "v")
-- Interpolate("psol2d", u, "p")


vanka = Vanka()
-- LineVanka(approxSpace)
-- vanka:set_num_steps(8,8,0,0,0,0)
-- vanka = DiagVanka()
vanka:set_damp(0.95)


vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-2, true))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-7)
baseConvCheck:set_reduction(1e-3)
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
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-2, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-2, true))


-- choose a solver
if discType=="stabil" then
--	solver = linSolver
else
--	solver = vankaSolver
end
solver = gmgSolver
solver = vankaSolver
-- solver = BiCGStabSolver

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

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element(VelCmp, "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")

step=0;
time=0;

-- create new grid function for old value
uOld = u:clone()

out:print("TimeDirichlet", u,0,0)

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
		
		-- compute CFL number 
		cflNumber(u,do_dt)

		print("++++++ TIMESTEP " .. step .. "  END ++++++");
		write("\n")
		l2error = L2Error("usol"..dim.."d", u, "u", time, 1, "Inner")
		write("L2Error in u component is "..l2error .."\n")
		l2error = L2Error("vsol"..dim.."d", u, "v", time, 1, "Inner")
		write("L2Error in v component is "..l2error .."\n")
		maxerror = MaxError("usol"..dim.."d", u, "u",time)
		write("Maximum error in u component is "..maxerror .."\n")
		maxerror = MaxError("vsol"..dim.."d", u, "v",time)
		write("Maximum error in v component is "..maxerror .."\n")
		write("\n")
		out:print("TimeDirichlet", u,step,time)
	end

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

print("done.")
