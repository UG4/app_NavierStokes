--------------------------------------------------------------------------------
--
--  Driven cavity problem
--
--  Author: Christian Wehner
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

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
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "COR", "Diffusion length type")
linred      = util.GetParam("-linred", 1e-1 , "Linear reduction")
nlintol     = util.GetParam("-nlintol", 1e-5, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")
boolStat = util.GetParamNumber("-stat", 1)
elemType = util.GetParam("-elem", "quad")
dt = util.GetParamNumber("-dt", 0.1)
reynoldsNr = util.GetParamNumber("-R",100)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "grids/unit_square_01_tri_unstruct_4bnd.ugx")
	else 
		gridName = util.GetParam("-grid", "grids/dc_quads.ugx") 
--		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2_four_bnd.ugx")	
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

print(" Chosen Parameters:")
print("    dim                 = " .. dim)
print("    numTotalRefs        = " .. numRefs)
print("    numPreRefs          = " .. numPreRefs)
print("    type                = " .. type)
if boolStat==false then
	print("    dt                  = " .. dt)
	print("    numTimeSteps        = " .. numTimeSteps)
end
print("    grid                = " .. gridName)
print("    v ansatz order      = " ..vorder)
print("    p ansatz order      = " ..porder)
print("    no laplace          = " .. tostring(bNoLaplace))
print("    exact jacobian      = " .. tostring(bExactJac))
print("    peclet blend        = " .. tostring(bPecletBlend))
print("    upwind              = " .. upwind)
print("    stab                = " .. stab)
print("    diffLength          = " .. diffLength)
print("    linear reduction    = " .. linred)
print("    linear tolerance    = " .. lintol)
print("    nonlinear reduction = " .. nlinred)
print("    nonlinear tolerance = " .. nlintol)
print("    stationary          = " .. boolStat)
print("    Reynolds-nr         = " .. reynoldsNr)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
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

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, type)

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(1.0/reynoldsNr);

--upwind if available
if type == "fv1" or type == "fvcr" then
	NavierStokesDisc:set_upwind(upwind)
	NavierStokesDisc:set_peclet_blend(bPecletBlend)
end

-- fv1 must be stablilized
if type == "fv1" then
	NavierStokesDisc:set_stabilization(stab, diffLength)
end

-- fe must be stabilized for (Pk, Pk) space
if type == "fe" and porder == vorder then
	NavierStokesDisc:set_stabilization(3)
end

--------------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------------

InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add({1,0}, "Top")

WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("Left,Right,Bottom")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)
domainDisc:add(WallDisc)

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

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
if type=="fvcr" then
	tBefore = os.clock()
	OrderCRCuthillMcKee(approxSpace,u,true)
--	CROrderCuthillMcKee(approxSpace,u,true,false,false,true)
	-- CROrderSloan(approxSpace,u,false,false,true)
	-- CROrderKing(approxSpace,u,true,false,false,true)
--	CROrderMinimumDegree(approxSpace,u,true,false,true)
	-- OrderLex(approxSpace, "lr");
	tAfter = os.clock()
	print("Ordering took " .. tAfter-tBefore .. " seconds.");
else
	if type=="fv1" then
--		OrderCuthillMcKee(approxSpace,true)
	end
end
u:set(0)

egsSolver = LinearSolver()
egsSolver:set_preconditioner(ElementGaussSeidel())
egsSolver:set_convergence_check(ConvCheck(10000, lintol, linred, true))

-- base solver
ilut = ILUT()
ilut:set_threshold(1e-7)
ilut:set_info(false)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(ilut)
ilutSolver:set_convergence_check(ConvCheck(100000,  lintol, linred, true))

if type=="fv1" then 
	ilutsmoother = ILUT()
	ilutsmoother:set_threshold(1e-5)
	smoother=ilutsmoother
elseif type=="fvcr" then 
	smoother=Vanka()
	smoother:set_damp(0.95)
elseif type=="fv" or type=="fe" then 
	smoother=ElementGaussSeidel()
end
smoother:set_damp(0.9)

if type=="fv1" then 
	basePre = ILUT()
	basePre:set_threshold(1e-7)
elseif type=="fvcr" then 
	basePre = Vanka()
elseif type=="fv" or type=="fe" then 
	basePre=ElementGaussSeidel()
end
baseSolver = LinearSolver()
baseSolver:set_preconditioner(basePre)
baseSolver:set_convergence_check(ConvCheck(10000, lintol*0.1,linred*0.1,false))

if type=="fv1" then   
--  for ilu 2 is sufficient
	numSmooth=2
else
--  for Vanka type smoothers there should be more smoothing steps on higher levels
	numSmooth=2*numRefs
end

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(baseSolver)
gmg:set_smoother(smoother)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(numSmooth)
gmg:set_num_postsmooth(numSmooth)
gmg:set_damp(MinimalResiduumDamping())
gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, lintol, linred, true))

BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(smoother)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, lintol, linred, true))

if type=="fv1" then 
	solver=gmgSolver
elseif type=="fvcr" then 
	solver=gmgSolver
elseif type=="fv" or type=="fe" then 
	solver=egsSolver
--	solver=ilutSolver
--	solver=LU()
--	solver=gmgSolver
--	solver=BiCGStabSolver
end

newtonSolver = NewtonSolver(op)
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(ConvCheck(10000, nlintol, nlinred, true))
newtonSolver:set_line_search(StandardLineSearch(30, 1.0, 0.9, true))
-- newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

SaveVectorForConnectionViewer(u, "StartSolution.vec")

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- start
time = 0.0
step = 0

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element(VelCmp, "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")

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
	
	out:print("DrivenCavity", u)
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
		out:print("TimeDrivenCavity", u)
	end
end
tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

dcevaluation(u,reynoldsNr)
