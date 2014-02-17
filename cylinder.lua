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
ug_load_script("util/domain_disc_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
numRefs 	= util.GetParamNumber("-numRefs", 0, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
bConvRates  = util.HasParamOption("-convRate", "compute convergence rates")

order 		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
vorder 		= util.GetParamNumber("-vorder", order, "order velocity space")
porder 		= util.GetParamNumber("-porder", order-1, "order pressure space")

bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type")
stab        = util.GetParam("-stab", "fields", "Stabilization type")
diffLength  = util.GetParam("-difflength", "raw", "Diffusion length type")

Viscosity	= util.GetParamNumber("-visco", 1e-3)
Um = 0.3

discType, vorder, porder = util.ns.parseParams()

if 	dim == 2 then gridName = util.GetParam("-grid", "grids/cylinder.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end


-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. discType)
print("    only stokes      = " .. tostring(bStokes))
print("    no laplace       = " .. tostring(bNoLaplace))
print("    exact jacobian   = " .. tostring(bExactJac))
print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)
print("    Viscosity	    = " .. Viscosity)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	
	local requiredSubsets = {}
	local dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
	
	return dom
end

function CreateApproxSpace(dom, discType, p)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, p, p-1)
	
	-- print statistic on the distributed dofs
	--approxSpace:init_levels()
	--approxSpace:init_top_surface()
	--approxSpace:print_statistic()
	--approxSpace:print_local_dof_statistic(2)
	
	return approxSpace
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

function CreateDomainDisc(approxSpace, discType, p)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		--NavierStokesDisc:set_pac_upwind(bPac)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		--NavierStokesDisc:set_stabilization(3)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(p*p+5)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(p*p+5)
	end
	
	-- setup Outlet
	--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	--OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		local H = 0.41
		return 4 * Um * y * (H-y) / (H*H), 0.0
	end
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
	--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	WallDisc:add("UpperWall,LowerWall,CylinderWall")
	
	-- Finally we create the discretization object which combines all the
	-- separate discretizations into one domain discretization.
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
	domainDisc:add(WallDisc)
	--domainDisc:add(OutletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

function CreateSolver(approxSpace, discType, p)

	local base = nil
	if discType == "fvcr" then
		base =  LinearSolver()
		base:set_preconditioner(DiagVanka())
		base:set_convergence_check(ConvCheck(10000, 1e-7, 1e-3, false))
	else
		base = SuperLU()
	end
	
	local smoother = nil
	if discType == "fvcr" then 
		smoother = Vanka()
	else
		local smooth = util.smooth.parseParams()
		smoother = util.smooth.create(smooth)
	end
	
	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	--gmg:set_damp(MinimalResiduumDamping())
	--gmg:set_damp(MinimalEnergyDamping())
	gmg:add_prolongation_post_process(AverageComponent("p"))
	--gmg:add_restriction_post_process(AverageComponent("p"))
	--gmg:set_debug(dbgWriter)
	
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))	
	end
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(ConvCheck(500, 1e-11, 1e-99, true))
	newtonSolver:set_line_search(StandardLineSearch(30, 1.0, 0.85, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)

	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
	AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

if bConvRates then
	
	local options = {	
	
		size = 				{12.5, 9.75}, -- the size of canvas (i.e. plot)
		sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
		font = 				"Arial",
		fontsize =			12,
		fontscale = 		1.4,
		
		logscale = 			true,
		grid = 				"lc rgb 'grey70' lt 0 lw 1", 
		linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
							linewidth = 3, pointsize = 1.3},
		border = 			" back lc rgb 'grey40' lw 2",
		decimalsign = 		",",
		key =	 			"on box lc rgb 'grey40' right bottom Left reverse spacing 1.5 width 1 samplen 2 height 0.5",
		tics =	 			{x = "nomirror out scale 0.75 format '%g' font ',8'",
							 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
		mtics =				5,
		slope = 			{dy = 3, quantum = 0.5, at = "last"},
		padrange = 			{ x = {0.8, 1.5}, y = {0.01, 1.5}},
	}

	if util.HasParamOption("-replot") then
		util.rates.static.replot(options); exit()
	end

	util.rates.static.compute(
	{	
		PlotCmps = { v = {"u","v"}, p = {"p"}},
		MeasLabel = function (disc, p) return disc.." $\\mathbb{Q}_{"..p.."}/\\mathbb{Q}_{"..(p-1).."}$" end,
		
		CreateDomain = CreateDomain,
		CreateApproxSpace = CreateApproxSpace,
		CreateDomainDisc = CreateDomainDisc,
		CreateSolver = CreateSolver,
		
		ComputeSolution = ComputeNonLinearSolution,
		
		DiscTypes = 
		{
		  {type = "fv", pmin = 3, pmax = 3, lmin = 1, lmax = numRefs},
		  --{type = "fe", pmin = 2, pmax = 5, lmin = 0, lmax = numRefs}
		},

		PrepareInitialGuess = function (u, lev, minLev, maxLev, domainDisc, solver)
			u[lev]:set(0.0)
		end,

		gpOptions = options,
		noplot = true,
		plotSol = true,
		MaxLevelPadding = function(p) return math.floor((p+1)/2) end,
		
	})
end

if not(bConvRates) then

	local p = vorder
	local dom = CreateDomain()
	local approxSpace = CreateApproxSpace(dom, discType, p)
	local domainDisc = CreateDomainDisc(approxSpace, discType, p)
	local solver = CreateSolver(approxSpace, discType, p)
	print(solver:config_string())
	
	local u = GridFunction(approxSpace)
	u:set(0)
	
	ComputeNonLinearSolution(u, domainDisc, solver)

	local D = 0.1
	local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, 3)
	print("Drag-Force: "..DL[1])
	print("Lift-Force: "..DL[2])
	print("C_D: "..(2*DL[1]/(1*math.pow(2*Um/3, 2)*D)))
	print("C_L: "..(2*DL[2]/(1*math.pow(2*Um/3, 2)*D)))

	PEval = GlobalGridFunctionNumberData(u, "p")
	Pa = PEval:evaluate({0.15, 0.2})
	Pe = PEval:evaluate( {0.25, 0.2} )

	print("Pa: "..Pa..", Pe: "..Pe..", DeltaP: "..Pa-Pe)

	local FctCmp = approxSpace:names()
	local VelCmp = {}
	for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
	vtkWriter = VTKOutput()
	vtkWriter:select(VelCmp, "velocity")
	vtkWriter:select("p", "pressure")
	vtkWriter:print("Cylinder", u)
end
