--------------------------------------------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--   Author: Josef Dubsky, Andreas Vogel
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim 		= util.GetParamNumber("-dim", 2)
type 		= util.GetParam("-type", "fvcr")
order 		= util.GetParamNumber("-order", 1)
vorder 		= util.GetParamNumber("-vorder", order)
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)
numRefs 	= util.GetParamNumber("-numRefs",2)
bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.GetParamNumber("-exactjac", 0)
R		 	= util.GetParamNumber("-R", 10)
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")


InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then gridName = util.GetParam("-grid", "grids/cylinder.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end


-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numTotalRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

requiredSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall", "CylinderWall"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"}; FctCmp = {"u", "p"};
elseif  dim == 2 then VelCmp = {"u", "v"}; FctCmp = {"u", "v", "p"};
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"};
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

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
elseif type == "fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else print("Disc Type '"..type.."' not supported.") exit() end

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(1.0/R)
NavierStokesDisc:set_peclet_blend(bPecletBlend)

if type == "fv1" then
	noUpwind = NavierStokesNoUpwind();
	fullUpwind = NavierStokesFullUpwind();
	skewedUpwind = NavierStokesSkewedUpwind();
	LPSUpwind = NavierStokesLinearProfileSkewedUpwind();
	POSUpwind = NavierStokesPositiveUpwind();
	
	fieldsStab = NavierStokesFIELDSStabilization()
	fieldsStab:set_upwind(fullUpwind)
	
	fieldsStab:set_diffusion_length("RAW")
	--fieldsStab:set_diffusion_length("FIVEPOINT")
	--fieldsStab:set_diffusion_length("COR")
	
	NavierStokesDisc:set_stabilization(fieldsStab)
end

if type == "fvcr" then
	noUpwind = NavierStokesNoUpwind();
	fullUpwind = NavierStokesFullUpwind();
	weightedUpwind = NavierStokesWeightedUpwind(0.6);
	NavierStokesDisc:set_upwind(weightedUpwind)
	
	NavierStokesDisc:set_upwind(weightedUpwind)
end

-- setup outlet
-- use "no normal stress" outflow condition
OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDisc:add("Outlet")

-- alternatively use "zero pressure" outflow condition
-- OutletDisc = DirichletBoundary()
-- OutletDisc:add(0.0, "p", "Outlet")

-- setup inlet
function inletVel2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return 4 * Um * y * (H-y) / (H*H), 0.0
end
InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Inlet")

-- setup wall
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("UpperWall,LowerWall,CylinderWall")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)
domainDisc:add(WallDisc)
domainDisc:add(OutletDisc)

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

u = GridFunction(approxSpace)

u:set(0)
time = 0.0

if type "fvcr" then
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
	gmg:set_base_linSolver(vankaBase)
	gmg:set_smoother(CRILUT)
	gmg:set_cycle_type(1)
	gmg:set_damp(MinimalResiduumDamping())
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	-- gmg:add_restriction_post_process(AverageComponent("p"))
	-- gmg:add_prolongation_post_process(AverageComponent("p"))
	
	--gmg:set_debug(dbgWriter)
	-- create Linear Solver
	BiCGStabSolver = BiCGStab()
	BiCGStabSolver:set_preconditioner(gmg)
	BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-1, true))
	
	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	gmgSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-1, true))
	
	-- choose a linSolver
	linSolver = gmgSolver
	-- linSolver = ilutSolver
	-- linSolver = vankaSolver
else
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_linSolver(LU())
	gmg:set_smoother(ILU())
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	
	-- create Linear Solver
	BiCGStabSolver = BiCGStab()
	BiCGStabSolver:set_preconditioner(gmg)
	BiCGStabSolver:set_convergence_check(ConvCheck(1000, 1e-7, 1e-3, true))
	
	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	gmgSolver:set_convergence_check(ConvCheck(1000, 1e-7, 1e-3, true))
	
	-- choose a linSolver
	linSolver = gmgSolver
end

-- Now we can set up the newton linSolver. We set the linear linSolver created above
-- as linSolver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver(domainDisc)
newtonSolver:set_linear_linSolver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(50, 1e-5, 1e-10, true))
newtonSolver:set_line_search(StandardLineSearch(20, 1.0, 0.5, true))
newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- Now we can apply the newton linSolver. A newton itertation is performed to find
-- the solution.
if newtonSolver:apply(u) == false then
	 print ("Newton linSolver apply failed."); exit();
end

-- Finally we're nearly done. The only thing left to do is to write the
-- solution to a file which can then be examined using e.g. Paraview.
-- (Open "Solution.vtu" in paraview to view the complete domain
vtkWriter = VTKOutput()
vtkWriter:select_all(false)
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "p")
vtkWriter:print("cylinder", u)