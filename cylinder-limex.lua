--------------------------------------------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--   Author: Arne Naegel
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

logAssistant = GetLogAssistant()
--logAssistant:set_debug_levels(10)
--logAssistant:set_debug_level("LIB_DISC_MULTIGRID", 10)
logAssistant:set_debug_level("LIB_LIMEX", 7)

local ARGS = {
  doSteadyState = util.HasParamOption("--steady-state", "Compute steady state solution"),
  bStokes   = util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed"),
  bExactJac   = util.HasParamOption("-exactjac", "If defined, exact jacobian used"),
  bLaplace  = not util.HasParamOption("-nolaplace", "If defined, only laplace term used"),
  bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used"),
  
  upwind      = util.GetParam("-upwind", "lps", "Upwind type"),
  
  stabGrad             = util.GetParamNumber("--stabGrad", 0.1, "Stabilization parameter."),
  stabStreamline      = util.GetParamNumber("--stabStreamline", 0.0, "Stabilization parameter."),
  stabDiv             = util.GetParamNumber("--stabDiv", 0.0, "Stabilization parameter."),

 
  limexDebugLevel = util.GetParamNumber("--limex-debug-level", 2, "debug level"),
  limexTOL = util.GetParamNumber("--limex-tol", 1e-2, "debug level"),
  limexNStages = util.GetParamNumber("--limex-num-stages", 4, "debug level"),
  
  solverID =  util.GetParam("--solver-id", "superlu", "superlu, gmg")
}

ARGS.doLimex = (ARGS.limexTOL > 0)

GetLogAssistant():set_debug_level("LIB_LIMEX", ARGS.limexDebugLevel)

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
numRefs 	= util.GetParamNumber("--numRefs", 1, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
bConvRates  = util.HasParamOption("-convuitRate", "compute convergence rates")
bBenchmarkRates = util.HasParamOption("-benchRate", "compute benchmark rates")

order 		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
vorder 		= util.GetParamNumber("-vorder", order, "order velocity space")
porder 		= util.GetParamNumber("-porder", order-1, "order pressure space")

stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type")

discType, vorder, porder = util.ns.parseParams()

local Viscosity	= 1e-3 --1e-3 --1e-3 -- 1.0 -- 1e-3 -- kinematic (nu=1/Re) or dynamic (mu) does not matter, since \rho=1. 
local Um = 0.3
if dim == 3 then Um = 0.45 end
local H = 0.41
local L = 0.1



-- Reference values for Schaefer /Turek benchmarks
local ref2D_1 = {
  CD = 5.57953523384,
  CL = 0.010618948146,
  DeltaP = 0.11752016697,
  
  Um =  1.5 --0.15 -- 1.5
}

local ref2D_3 = {
  
  CD = 2.950921575, tCD=3.93625, -- maximum value and time (cited according to in John, Rang)
  CL = 0.47795, tCL= 5.693125,
  DeltaP = -0.1116,  -- at t=8
  
  Um = 1.5
}


ref = ref2D_1
Um = ref.Um
local Umean2 = math.pow(2/3*Um, 2)


if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/cylinderg.ugx")
	--gridName = util.GetParam("-grid", "grids/box.ugx")
	--gridName = util.GetParam("-grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_rotate_box_tri_fine.ugx")
elseif dim == 3 then
	gridName = util.GetParam("-grid", "grids/cylinder3d.ugx")
--	gridName = util.GetParam("-grid", "grids/cylinder3d_fine.ugx")
else print("Selected Dimension not supported. Exiting."); exit(); end


-- Lets write some info about the choosen parameter
print(" Selected Parameter:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. discType)
print("    only stokes      = " .. tostring(ARGS.bStokes))
print("    only laplace       = " .. tostring(ARGS.bLaplace))
print("    exact jacobian   = " .. tostring(ARGS.bExactJac))
print("    peclet blend     = " .. tostring(ARGS.bPecletBlend))
print("    upwind           = " .. ARGS.upwind)
print("    stab             = " .. stab)
print("    stabGrad         = " .. ARGS.stabGrad)
print("    stabDiv          = " .. ARGS.stabDiv)
print("    stabStreamline   = " .. ARGS.stabStreamline)
print("    diffLength       = " .. diffLength)
print("    LIMEX.TOL       = " .. ARGS.limexTOL)
print("    LIMEX.NumStages       = " .. ARGS.limexNStages)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	local dom = Domain()
	LoadDomain(dom, gridName)
	
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner =  GlobalDomainRefiner(dom)
	
	write("Pre-Refining("..numPreRefs.."): ")
	for i=1,numPreRefs do write(i .. " ");	refiner:refine(); end
	write("done. Distributing...")
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end
	write("done.\n")
	
	SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
	
	return dom
end

function CreateApproxSpace(dom, discType, vorder, porder)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder)
	
	
	-- print statistic on the distributed dofs
  approxSpace:init_levels()
  approxSpace:init_top_surface()

	

	
	approxSpace:print_statistic()
	approxSpace:print_local_dof_statistic(2)
	
	
--  OrderLex(approxSpace,"x")
	-- OrderCuthillMcKee(approxSpace,true)
	return approxSpace
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------
globalNSDisc = nil

Cylinder2D_FE = {

  disc = {
    walls = "UpperWall,LowerWall,CylinderWall",
    inlet = function (x, y, t) return 4 * Um * y * (H-y) / (H*H), 0.0 end,
  }

}

Cylinder3D_FE = {
 
  disc = {
    walls = "UpperWall,LowerWall,CylinderWall,FrontWall,BackWall",
    inlet = function (x, y, z, t) return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0  end,
  }
  
}


-- Creates the domain discretization
function CreateDomainDisc(approxSpace, discType) --, vorder, porder)

	local FctCmp = approxSpace:names()
	local NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(ARGS.bExactJac)
-- 	NavierStokesDisc:set_exact_jacobian(1.0)
	NavierStokesDisc:set_stokes(ARGS.bStokes)
	NavierStokesDisc:set_laplace(ARGS.bLaplace)
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	globalNSDisc = NavierStokesDisc
				
	local porder = approxSpace:lfeid(dim):order()
	local vorder = approxSpace:lfeid(0):order()
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(ARGS.upwind)
		NavierStokesDisc:set_peclet_blend(ARGS.bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(true)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if (discType == "fe") and (porder == vorder) then
		NavierStokesDisc:set_stabilization(ARGS.stabGrad)
		NavierStokesDisc:set_stab_streamline(ARGS.stabStreamline) -- 1.0
		
	end
	if discType == "fe" then
	 NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	 -- NavierStokesDisc:set_quad_order(3)
	 NavierStokesDisc:set_stab_div(ARGS.stabDiv)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	
	-- setup Outlet
	-- OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	-- OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		return 4 * Um * y * (H-y) / (H*H)*math.sin(math.pi*t/8.0), 0.0
	end
	function inletVelX2d(x, y, t)
    return 4 * Um * y * (H-y) / (H*H)*math.sin(math.pi*t/8.0)
  end
	function inletVel3d(x, y, z, t)
		return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0
	end
	
	local InletDisc = NavierStokesInflow(NavierStokesDisc)
	--InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet") -- 
	InletDisc:add("inletVel"..dim.."d", "Inlet")
	
	-- John's (physically unrealistic) BC
  --local InletDisc = DirichletBoundary()
	--InletDisc:add("inletVelX"..dim.."d", "u" ,"Inlet, Outlet")
	--InletDisc:add(0.0, "v" ,"Inlet, Outlet")
	
	--setup Walls
	local WallDisc = NavierStokesWall(NavierStokesDisc)
	if dim == 2 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall")
	elseif dim == 3 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall")	
	end

	local DirichletBnd = DirichletBoundary()
  -- DirichletBnd:add(0, "p", "FIXP")  -- fix pressure => oscillations
  -- DirichletBnd:add("inletVelX2d", "u", "Outlet")  -- fix pressure
  -- DirichletBnd:add(0.0, "v", "Outlet")  -- fix pressure
	
	local NoSlipBnd = DirichletBoundary()
	NoSlipBnd:add(0, "u", "UpperWall,LowerWall,CylinderWall")  -- no slip
	NoSlipBnd:add(0, "v", "UpperWall,LowerWall,CylinderWall")  -- no slip
	
	
	-- Finally we create the discretization object which combines all the
	-- separate discretizations into one domain discretization.
	local domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
	--domainDisc:add(WallDisc)
	domainDisc:add(DirichletBnd)
	domainDisc:add(NoSlipBnd)
	--domainDisc:add(OutletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------
-- 
function CreateSolver(approxSpace, discType, p)

	local base = SuperLU()
	
	local smoother = nil
	if discType == "fvcr" or discType == "fecr" then 
		smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	elseif discType == "fv1" then 
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		 smoother = ComponentGaussSeidel(1.0, {"p"})
		 smoother:set_alpha(1.0)
		 smoother:set_beta(1.0)
		 smoother:set_weights(true)
	end
	
	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local cycleType = "W"
	local bRAP= true
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycleType, base, baseLev, bRAP)
	--gmg:add_prolongation_post_process(AverageComponent("p"))
	local transfer = StdTransfer()
	transfer:enable_p1_lagrange_optimization(false)
	gmg:set_transfer(transfer)
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	if ARGS.bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))	
	end
		
	local convCheck = ConvCheck(50, 1e-11, 1e-99, true)
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
	-- newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end




function ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)

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
		
		ComputeSolution = ComputeNonLinearSteadyStateSolution,
		
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

	local dom = CreateDomain()
	local plots = {}	

	 function ComputeSpace(discType, p, ppress, minLev, maxLev)

		local file = table.concat({"dc",discType,p,ppress},"_")..".dat"
		--[[
		local discLabel = discType.." Q_"..p.."/Q_"..(ppress)
		local CDLabel = "|C_D - C_D_h|"
		local CLLabel = "|C_L - C_L_h|"
		local DeltaPLabel = "|P - P_L_h|"
		--]]
		local discLabel = discType.." $\\mathbb{Q}_{"..p.."}/\\mathbb{Q}_{"..(ppress).."}$"
		local CDLabel = "$|c_{w,h} - c_w^{\\text{Ref}}|$"
		local CLLabel = "$|c_{a,h} - c_a^{\\text{Ref}}|$"
		local DeltaPLabel = "$|\\Delta_{p,h} - \\Delta_p^{\\text{Ref}}|$"
		
		local function addPlot(name, dataset, label)
			plots[name] = plots[name] or {}
			table.insert( plots[name], dataset)			
			plots[name].label = label			
		end
		
		addPlot("CD_DoF", {label=discLabel, file=file, style="linespoints", 1, 3},
				{ x = "Anzahl Unbekannte", y = CDLabel})

		addPlot("CL_DoF", {label=discLabel, file=file, style="linespoints", 1, 4},
				{ x = "Anzahl Unbekannte", y = CLLabel})

		addPlot("DeltaP_DoF", {label=discLabel, file=file, style="linespoints", 1, 5},
				{ x = "Anzahl Unbekannte", y = DeltaPLabel})

		addPlot("CD_h", {label=discLabel, file=file, style="linespoints", 2, 3},
				{ x = "h (Gitterweite)", y = CDLabel})

		addPlot("CL_h", {label=discLabel, file=file, style="linespoints", 2, 4},
				{ x = "h (Gitterweite)", y = CLLabel})

		addPlot("DeltaP_h", {label=discLabel, file=file, style="linespoints", 2, 5},
				{ x = "h (Gitterweite)", y = DeltaPLabel})

		if not util.HasParamOption("-replot") then

			local approxSpace = util.ns.CreateApproxSpace(dom, discType, p, ppress)
			local domainDisc = CreateDomainDisc(approxSpace, discType, p, ppress)
			local solver = CreateSolver(approxSpace, discType, p)
		
			local h, DoFs, level = {}, {}, {}	
			local meas = {CD = {}, CL = {}, DeltaP = {}}
			for n, _ in pairs(meas) do
				meas[n] = {value = {}, prev = {error = {}, rate = {}}, exact = {error = {}, rate = {}}}
			end	
			
			local uPrev = nil
			for lev = minLev, maxLev do
				write("\n>> Computing Level "..lev..", "..discType..", "..p..", "..ppress..".\n")
		
				local u = GridFunction(approxSpace, lev)

				if uPrev ~= nil and discType ~= "fvcr" and discType ~= "fecr" then	
					Prolongate(u, uPrev);
					AdjustMeanValue(u, "p")
					u:enforce_consistent_type()
					u:check_storage_type()
				else
					u:set(0)	
				end
				
				if lev > minLev + 2 then globalNSDisc:set_exact_jacobian(true) end
				
				ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)
				u:check_storage_type()

				local FctCmp = approxSpace:names()
				local VelCmp = {}
				
				for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
				local vtkWriter = VTKOutput()
				vtkWriter:select(VelCmp, "velocity")
				vtkWriter:select("u", "u")
				vtkWriter:select("v", "v")
				vtkWriter:select("p", "pressure")
				vtkWriter:print("Cylinder"..table.concat({discType,p,ppress,"lev",lev},"_"), u)
		
				-- h/DoF statistic
				DoFs[lev] = u:num_dofs()
				h[lev] =  MaxElementDiameter(dom, lev) 
				level[lev] = lev		
				
				-- C_D / C_L
				local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p*p+5)
				meas.CD.value[lev] = 2*DL[1]/(Umean2*L)
				meas.CL.value[lev] = 2*DL[2]/(Umean2*L)
			
				-- Delta_P
				local PEval = GlobalGridFunctionNumberData(u, "p")
				meas.DeltaP.value[lev] = PEval:evaluate(Vec2d(0.15, 0.2)) - PEval:evaluate(Vec2d(0.25, 0.2))
		
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
		
				write("\n>> Stats: lev "..lev..", "..discType..", "..p..", "..ppress..".\n")
				table.print(values, {heading = heading, format = format, 
									 hline = true, vline = true, forNil = "--"})
									 
				uPrev = u
				--collectgarbage()
			end
			
			local cols = {DoFs, h, meas.CD.exact.error, meas.CL.exact.error, meas.DeltaP.exact.error}
			gnuplot.write_data(file, cols)	
			
			approxSpace, domainDisc, solver = nil, nil, nil
			collectgarbage()
		end
	end
	
--	ComputeSpace("fecr", 1, 0, numPreRefs, numRefs)	
--	ComputeSpace("fvcr", 1, 0, numPreRefs, numRefs)	
--	ComputeSpace("fe", 1, 1, numPreRefs, numRefs)	
--	ComputeSpace("fe", 2, 1, numPreRefs, numRefs-1)	
--	ComputeSpace("fe", 3, 2, numPreRefs, numRefs-2)	
	ComputeSpace("fv1", 1, 1, numPreRefs, numRefs)	
--	ComputeSpace("fv", 2, 1, numPreRefs, numRefs-1)	
--	ComputeSpace("fv", 3, 2, numPreRefs, numRefs-2)	
	
	if util.HasParamOption("-replot") then
		local texOptions = {	
		
			size = 				{12.5, 6.75}, -- the size of canvas (i.e. plot)
			sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
			font = 				"Arial",
			fontsize =			12,
			fontscale = 		0.7,
			
			logscale = 			true,
			grid = 				"lc rgb 'grey70' lt 0 lw 1", 
			linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
								linewidth = 3, pointsize = 1.3},
			border = 			" back lc rgb 'grey40' lw 2",
			decimalsign = 		",",
			key =	 			"on box lc rgb 'grey40' left bottom Left reverse spacing 2 width 1.1 samplen 2 height 0.5",
			tics =	 			{x = "nomirror out scale 0.75 format '%.te%01T' font ',8'",
								 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
			mtics =				5,
			slope = 			{dy = 3, quantum = 0.5, at = "last"},
			padrange = 			{ x = {0.6, 2}, y = {0.6, 2}},
		}

		local pdfOptions = {	
		
			size = 				{12.5, 9.75}, -- the size of canvas (i.e. plot)
			sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
			font = 				"Arial",
			fontsize =			8,
			
			logscale = 			true,
			grid = 				"lc rgb 'grey70' lt 0 lw 1", 
			linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
								linewidth = 3, pointsize = 1.3},
			border = 			" back lc rgb 'grey40' lw 2",
			decimalsign = 		",",
			key =	 			"on box lc rgb 'grey40' left bottom Left reverse spacing 2 width 1.1 samplen 2 height 0.5",
			tics =	 			{x = "nomirror out scale 0.75 format '%g' font ',8'",
								 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
			mtics =				5,
			slope = 			{dy = 3, quantum = 0.25, at = "last"},
			padrange = 			{ x = {0.6, 10}, y = {0.1, 1.1}},
		}
	
		for name,data in pairs(plots) do
			gnuplot.plot(name..".pdf", data, pdfOptions)
			gnuplot.plot(name..".tex", data, texOptions)
		end
	end
end

if not(bConvRates) and not(bBenchmarkRates) then

	local p = vorder
	local dom = CreateDomain()
	local approxSpace = CreateApproxSpace(dom, discType, p, porder)
	local domainDisc = CreateDomainDisc(approxSpace, discType, p)
	local solver = CreateSolver(approxSpace, discType, p)

  local dbgWriter = GridFunctionDebugWriter(approxSpace)
  -- dbgWriter.set_conn_viewer_output(true)

 -- local transfer = 
  local transfer = StdTransfer()
  transfer:enable_p1_lagrange_optimization(true)
  -- transfer:set_debug(dbgWriter)
  
  
	local solverDesc = {
  
  type = "newton",
  lineSearch = StandardLineSearch(10, 1.0, 0.5, true, true),
  convCheck = ConvCheck(50, 1e-13, 1e-8, true, false),-- "standard",  -- verbose, suppress
  
  
   linSolver = {
    type = "bicgstab", --"bicgstab",
    -- type = "superlu", --"bicgstab",
    convCheck = ConvCheck(100, 1e-16, 1e-8, true),
   }
  }
 
  if (ARGS.solverID == "superlu") then
    solverDesc.linSolver.type = "superlu"
   elseif (ARGS.solverID == "gmg") then
   
    solverDesc.linSolver.type = "linear"
    solverDesc.linSolver.precond =  {
      type    = "gmg",
      
      approxSpace = approxSpace,
      baseLevel = 0,
      baseSolver  = "superlu",
      cycle = "W",
      
      rap = true,
      -- discretization = domainDisc,
      
      preSmooth = 1,
      postSmooth = 1,

      smoother  =  -- [[jac, egs,
      {
        type = "egs",
        damping = 1.0,
        vertex = {{"p"}, {"u", "v"}},
 
      },--]]
      
      
      -- transfer = transfer, 
     debug = true,
    
      -- transfer = "std",
      
      -- [[ 
      debugSolver = {
      
        type = "linear",
        convCheck = ConvCheck(20, 5e-12, 1e-99, true),
        precond =  {
          type = "ssc",
          vertex = {{"p"}, {"u","v"}}
        },
        --[[egs,
       {
            type = "egs",
            damping = 1.0,
            vertex = {{"p"}, {"u", "v"}}
     
        },--]]
        
        approxSpace = approxSpace,
        -- debug = true
       }, -- debugSolver
      --]]
     }  
  
 
 else  


     precond = {
        type = "ssc",
        damping = 1.0,
        vertex = {{"p"}, {"u", "v"}}
     } --]]
 
    end
 
 
 
 
  
 -- util.debug_writer = dbgWriter

	local solver = util.solver.CreateSolver(solverDesc, solverutil)
	 --solver:set_debug(GridFunctionDebugWriter(approxSpace))
  print(solver:config_string())




  
	 -- VTK writer.
  local FctCmp = approxSpace:names()
  local VelCmp = {}
  for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
  
  local vtkWriter = VTKOutput()
  vtkWriter:select(VelCmp, "velocity")
  vtkWriter:select("p", "pressure")
  
  
	-- Create grid function.
	local u = GridFunction(approxSpace)
	u:set(0)
	
	if (ARGS.doSteadyState) then
	   -- Steady state solution.
	   ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)
	   vtkWriter:print("CylinderSteadyState", u)
  else
    
      
      -- Transient solution.
      local cTransient = 
      {
        -- Start and stop time --
        tStart = 0.0,
        tStop = 8.0,

        -- Fractional-step-theta (same amount of work as CN, but stable) --
        scheme = "fracstep", maxStepSize = 0.1, minStepSize = 8.0/4096.0, redStepSize = 0.5,

      --  scheme = "limex", maxStepSize = 0.5, minStepSize = (8.0/2048.0)*1e-4, redStepSize = 0.5, 
       

        -- Implicit Euler (requires a small time step) --
        --scheme = "impleuler", maxStepSize = 0.001, minStepSize = 0.0001, redStepSize = 0.5

        -- SDIRK (3rd-order) a.k.a. Alexander3 (we call it `sdirk3`) --
        --scheme = "sdirk", orderOrTheta = 3, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 

        -- SDIRK (4th-order) a.k.a. Hairer, Wanner, L-stable DIRK (we call it `sdirk4`), TODO: Check --
        --scheme = "sdirk", orderOrTheta = 4, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 
    }
  
    local ComputeEffectiveQuantities = nil
    
    if ARGS.doLimex then 
      cTransient.scheme = "limex" 
      cTransient.maxStepSize = 0.5
      cTransient.minStepSize = (8.0/2048.0)*1e-8
      cTransient.redStepSize = 0.5
    end
    
    
-- Evaluate drag, lift and deltaP
function EvalIntegralQuantities2D (u, step, time)
    print("EvalIntegralQuantities2D")

    local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
    local C_D = 2*DL[1]/(Umean2*L)
    local C_L = 2*DL[2]/(Umean2*L)
  
    local PEval = GlobalGridFunctionNumberData(u, "p")
    local Delta_P = PEval:evaluate(Vec2d(0.15, 0.2)) - PEval:evaluate(Vec2d(0.25, 0.2))
  
    print("EVAL_P1:\t"..time.."\t"..PEval:evaluate(Vec2d(0.15, 0.2)))
    print("EVAL_P2:\t"..time.."\t"..PEval:evaluate(Vec2d(0.25, 0.2)))
    print("EVAL_DELTA_P:\t"..time.."\t"..Delta_P)
  
    print("EVAL_C_D:\t"..time.."\t"..C_D)
    print("EVAL_C_L:\t"..time.."\t"..C_L)
end
    
    
    ComputeEffectiveQuantities = EvalIntegralQuantities2D
   
    if (cTransient.scheme == "limex") then 
    
      RequiredPlugins({"Limex"})
      ug_load_script("plugins/Limex/limex_util.lua")
    
      local adaptiveStepInfo = {}
      adaptiveStepInfo["STAGES"] = ARGS.limexNStages -- tolerance
      adaptiveStepInfo["TOLERANCE"] = ARGS.limexTOL-- tolerance
      adaptiveStepInfo["REDUCTION"] = 0.5  -- reduction
      adaptiveStepInfo["INCREASE"]  = 2.0 -- increase of time step
      adaptiveStepInfo["SAFETY"] = 0.8    -- safety factor
  -- local errorEst = adaptiveStepInfo["ESTIMATOR"]
       
       adaptiveStepInfo["SPACES"] = {
        --  L2ComponentSpace("p", 2), 
         H1SemiComponentSpace("u", 4), 
         H1SemiComponentSpace("v", 4), -- Viscosity
        --  L2QuotientSpace("p", 2), -- need to factor out constants!
       }
       
       -- adaptiveStepInfo["DEBUG"] = GridFunctionDebugWriter(approxSpace) -- enable output of matrices
       
    util.SolveNonlinearProblemLimex(
        u, domainDisc, solver,
        vtkWriter, "limex_solution",
        cTransient.tStart, cTransient.tStop, 
        cTransient.maxStepSize*0.01,  -- step size
        cTransient.minStepSize, cTransient.maxStepSize, 
        adaptiveStepInfo, ComputeEffectiveQuantities)
    else
    
     util.SolveNonlinearTimeProblem(
        u, domainDisc, solver,
        vtkWriter, "std_solution",
        cTransient.scheme, 
        cTransient.orderOrTheta, 
        cTransient.tStart, cTransient.tStop, 
        cTransient.maxStepSize, cTransient.minStepSize, cTransient.redStepSize, 
        false, false, ComputeEffectiveQuantities)
    
    end
    
        
	end
	


  -- Compute drag and lift coefficients.
	if dim == 2 then 
		EvalIntegralQuantities2D(u)
	end	
end
