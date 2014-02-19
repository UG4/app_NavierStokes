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
bBenchmarkRates = util.HasParamOption("-benchRate", "compute benchmark rates")

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

discType, vorder, porder = util.ns.parseParams()

local Viscosity	= 1e-3
local Um = 0.3
local H = 0.41
local L = 0.1
local Umean2 = math.pow(2/3*Um, 2)

local ref = {}
ref.CD = 5.57953523384
ref.CL = 0.010618948146
ref.DeltaP = 0.11752016697

if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/cylinder.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder-rims.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder-rims2.ugx")
	gridName = util.GetParam("-grid", "grids/box.ugx")
	--gridName = util.GetParam("-grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_rotate_box_tri_fine.ugx")
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

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	
	-- create Instance of a Domain
	local dom = Domain()
	
	-- load domain
	write("Loading Domain "..gridName.." ... ") 
	LoadDomain(dom, gridName)
	write("done. ")
	
	-- create Refiner
	if numPreRefs > numRefs then
		print("numPreRefs must be smaller than numRefs. Aborting.");
		exit();
	end
	
	if numPreRefs > numRefs then
		numPreRefs = numRefs
	end
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner =  GlobalDomainRefiner(dom)
	local refProjector = DomainRefinementProjectionHandler(dom)
	--refProjector:set_callback("CylinderWall", SphereProjector(dom, 0.2, 0.2, 0, 0.05))
	--refProjector:set_callback("CylinderRim1", SphereProjector(dom, 0.2, 0.2, 0, 0.0525))
	--refProjector:set_callback("CylinderRim2", SphereProjector(dom, 0.2, 0.2, 0, 0.05625))
	--refProjector:set_callback("CylinderRim3", SphereProjector(dom, 0.2, 0.2, 0, 0.063367748))
	--refProjector:set_callback("Inner", SubdivisionLoopProjector(dom))
	refiner:set_refinement_callback(refProjector)
	
	write("Pre-Refining("..numPreRefs.."): ")
	-- Performing pre-refines
	for i=1,numPreRefs do
		write(i .. " ")
		refiner:refine()
	end
	write("done. Distributing...")
	-- Distribute the domain to all involved processes
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	
	-- Perform post-refine
	for i=numPreRefs+1,numRefs do
		refiner:refine()
		write(i-numPreRefs .. " ")
	end
	write("done.\n")
	
	-- Now we loop all subsets an search for it in the SubsetHandler of the domain
	if neededSubsets ~= nil then
		if util.CheckSubsets(dom, neededSubsets) == false then 
			print("Something wrong with required subsets. Aborting.");
			exit();
		end
	end
	
	
	--clean up
	if refiner ~= nil then
		delete(refiner)
	end
	
	-- return the created domain
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
	transfer = StdTransfer()
	transfer:enable_p1_lagrange_optimization(false)
	gmg:set_transfer(transfer)
	
--	LinearIteratorProduct({smoother, gmg})
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, LinearIteratorProduct({gmg, smoother}))
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 5e-3, true))	
	end
		
	local convCheck = ConvCheck(500, 1e-11, 1e-99, true)
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(5, 1.0, 0.9, true, true))
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
		  --{type = "fv", pmin = 3, pmax = 3, lmin = 1, lmax = numRefs},
		  {type = "fe", pmin = 2, pmax = 2, lmin = 0, lmax = numRefs}
		},

		PrepareInitialGuess = function (u, lev, minLev, maxLev, domainDisc, solver)
			u[lev]:set(0.0)
		end,

		gpOptions = options,
		noplot = true,
		plotSol = true,
		--MaxLevelPadding = function(p) return math.floor((p+1)/2) end,
		MaxLevelPadding = function(p) return 0 end,
		
	})
end

if bBenchmarkRates then

	local p = vorder
	local dom = CreateDomain()
	local approxSpace = CreateApproxSpace(dom, discType, p)
	local domainDisc = CreateDomainDisc(approxSpace, discType, p)
	local solver = CreateSolver(approxSpace, discType, p)

	local h, DoFs, level = {}, {}, {}	
	local meas = {CD = {}, CL = {}, DeltaP = {}}
	for n, _ in pairs(meas) do
		meas[n] = {value = {}, prev = {error = {}, rate = {}}, exact = {error = {}, rate = {}}}
	end	
	
	local uPrev = nil
	local minLev = 0
	local maxLev = numRefs
	for lev = minLev, maxLev do
		write("\n>> Computing Level "..lev..", "..discType..", "..p..".\n")

		local u = GridFunction(approxSpace, lev)
		
		if uPrev ~= nil then	
			Prolongate(u, uPrev);
			AdjustMeanValue(u, "p")
		else
			u:set(0)	
		end
		u:set(0)	
		
		write(">> Start: Computing solution on level "..lev..".\n")
		local convCheck = CompositeConvCheck(approxSpace, 500, 1e-11, 1e-99)
		convCheck:set_all_component_check(1e-11, 1e-99)
		convCheck:set_level(lev)
		solver:set_convergence_check(convCheck)
		--solver:set_debug(GridFunctionDebugWriter(approxSpace))
			
		ComputeNonLinearSolution(u, domainDisc, solver)
		write(">> End: Solver done.\n")

		-- h/DoF statistic
		DoFs[lev] = u:size()
		h[lev] =  MaxElementDiameter(dom, lev) 
		level[lev] = lev		
		
		-- C_D / C_L
		local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
		meas.CD.value[lev] = 2*DL[1]/(Umean2*L)
		meas.CL.value[lev] = 2*DL[2]/(Umean2*L)
	
		-- Delta_P
		local PEval = GlobalGridFunctionNumberData(u, "p")
		meas.DeltaP.value[lev] = PEval:evaluate({0.15, 0.2}) - PEval:evaluate( {0.25, 0.2} )

		for n, _ in pairs(meas) do
			local quant = meas[n]
			quant.exact.error[lev] = math.abs(quant.value[lev] - ref[n])
			
			if lev > minLev then
				quant.prev.error[lev-1] = math.abs(quant.value[lev] - quant.value[lev-1])
			end			
			for _, t in ipairs({"exact", "prev"}) do
				local type = quant[t]
				if type.error[lev-2] ~= nil and type.error[lev-1] ~= nil then
					local fac = type.error[lev-2] / type.error[lev-1]
					type.rate[lev-1] = math.log(fac) / math.log(2) --math.log(h[lev-1]/h[lev])
				end			
				if type.error[lev-1] ~= nil and type.error[lev] ~= nil then
					local fac = type.error[lev-1] / type.error[lev]
					type.rate[lev] = math.log(fac) / math.log(2) --math.log(h[lev-1]/h[lev])
				end			
			end
		end
				
		-- print
		local values = {level, h, DoFs}
		local heading = {"L", "h", "#DoFs"}
		local format = {"%d", "%.2e", "%d"}

		for n, _ in pairs(meas) do
			local quant = meas[n]
			table.append(values, {quant.value}) 
			table.append(heading,{n})
			table.append(format, {"%.8f"})

			for _, t in ipairs({"exact", "prev"}) do
				local type = quant[t]
				table.append(values, {type.error, type.rate}) 
				table.append(heading,{n.." "..t, "rate"})
				table.append(format, {"%.3e", "%.3f"})
			end
		end

		table.print(values, {heading = heading, format = format, 
							 hline = true, vline = true, forNil = "--"})
							 
		 uPrev = u
		collectgarbage()
	end
	
	approxSpace, domainDisc, solver = nil, nil, nil
	collectgarbage()
	
end

if not(bConvRates) and not(bBenchmarkRates) then

	local p = vorder
	local dom = CreateDomain()
	local approxSpace = CreateApproxSpace(dom, discType, p)
	local domainDisc = CreateDomainDisc(approxSpace, discType, p)
	local solver = CreateSolver(approxSpace, discType, p)
	--solver:set_debug(GridFunctionDebugWriter(approxSpace))
			
	print(solver:config_string())
	
	local u = GridFunction(approxSpace)
	u:set(0)
	
	ComputeNonLinearSolution(u, domainDisc, solver)

	local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
	local C_D = 2*DL[1]/(Umean2*L)
	local C_L = 2*DL[2]/(Umean2*L)

	local PEval = GlobalGridFunctionNumberData(u, "p")
	local Delta_P = PEval:evaluate({0.15, 0.2}) - PEval:evaluate( {0.25, 0.2} )

	print("p1: "..PEval:evaluate({0.15, 0.2}))
	print("p2: "..PEval:evaluate({0.25, 0.2}))

	print("C_D - ref.CD: "..string.format("%.3e", C_D - ref.CD))
	print("C_L - ref.CL: "..string.format("%.3e", C_L - ref.CL))
	print("Delta_P - ref.DeltaP: "..string.format("%.3e", Delta_P - ref.DeltaP))
	

	local FctCmp = approxSpace:names()
	local VelCmp = {}
	for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
	vtkWriter = VTKOutput()
	vtkWriter:select(VelCmp, "velocity")
	vtkWriter:select("p", "pressure")
	vtkWriter:print("Cylinder", u)
end
