--------------------------------------------------------------------------------
--
--  Driven cavity problem
--
--  Author: Christian Wehner
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "fvcr" or "fv1"
discType = util.GetParam("-type", "fvcr")
boolStat = util.GetParamNumber("-stat", 1)
elemType = util.GetParam("-elem", "quad")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then gridName = util.GetParam("-grid", "grids/unit_square_01_tri_unstruct_4bnd.ugx")
	else 
		gridName = util.GetParam("-grid", "grids/dc_quads.ugx") 
--		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2_four_bnd.ugx")	
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

dt = util.GetParamNumber("-dt", 0.1)
reynoldsNr = util.GetParamNumber("-Re",5000)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

print(" Chosen Parameters:")
print("    dim          = " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type         = " .. discType)
if boolStat==false then
print("    dt           = " .. dt)
print("    numTimeSteps = " .. numTimeSteps)
end
print("    grid      	= " .. gridName)
print("    stationary   = " .. boolStat)
print("    Reynolds-nr  = " .. reynoldsNr)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"}; FctCmp = {"u", "p"};
elseif  dim == 2 then VelCmp = {"u", "v"}; FctCmp = {"u", "v", "p"};
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"};
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

if discType=="fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else
	approxSpace:add_fct(FctCmp, "Lagrange", 1)
end

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)

if discType=="fvcr" then
	noUpwind = NavierStokesNoUpwind();
	fullUpwind = NavierStokesFullUpwind();
	weightedUpwind = NavierStokesWeightedUpwind(0.5);
	NavierStokesDisc:set_conv_upwind(fullUpwind)
else	
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
NavierStokesDisc:set_defect_upwind(true)
NavierStokesDisc:set_kinematic_viscosity(1.0/reynoldsNr);
NavierStokesDisc:set_source({0,0})

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
if discType=="fvcr" then
	tBefore = os.clock()
	-- OrderCRCuthillMcKee(approxSpace,u,true)
	CROrderCuthillMcKee(approxSpace,u,true,false,false,true)
	-- CROrderSloan(approxSpace,u,false,false,true)
	-- CROrderKing(approxSpace,u,true,false,false,true)
--	CROrderMinimumDegree(approxSpace,u,true,false,true)
	-- OrderLex(approxSpace, "lr");
	tAfter = os.clock()
	print("Ordering took " .. tAfter-tBefore .. " seconds.");
else
	OrderCuthillMcKee(approxSpace,true)
end
u:set(0)

vanka = Vanka()
vanka:set_damp(0.95)

crilut = CRILUT()
crilut:set_threshold(1e-1,1e-3,1e-3,1e-4)
-- crilut:set_threshold(1e-4)
crilut:set_info(true)
crilut:set_damp(1)

crilutSolver = LinearSolver()
crilutSolver:set_preconditioner(crilut)
crilutSolver:set_convergence_check(ConvCheck(10000, 1e-5, 2e-1, true))

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(100000, 1e-5, 2e-1, true))

vankaBase = LinearSolver()
vankaBase:set_preconditioner(Vanka())
vankaBase:set_convergence_check(ConvCheck(10000, 1e-5, 1e-2, false))

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(vankaBase)
gmg:set_smoother(vanka)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(5)
gmg:set_num_postsmooth(5)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(0.8)
-- gmg:set_damp(MinimalEnergyDamping())

if discType=="fv1" then
	crilut:set_threshold(1e-5)
	gmg:set_base_solver(crilutSolver)
	gmg:set_smoother(crilut)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
end

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-5, 2e-1, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-5, 2e-1, true))

-- choose a solver
solver = BiCGStabSolver
solver = gmgSolver
solver = crilutSolver
-- solver = vankaSolver

newtonSolver = NewtonSolver(op)
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(ConvCheck(10000, 1e-5, 1e-10, true))
newtonSolver:set_line_search(StandardLineSearch(20, 1.0, 0.9, true))
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
	
	out:print("DCSolution", u)
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
		out:print("timeDCSolution", u)
	end
end
tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

dcevaluation(u,reynoldsNr)
