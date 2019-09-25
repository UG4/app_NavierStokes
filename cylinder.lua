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
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type")

local discType, vorder, porder = util.ns.parseParams()

local Viscosity	= 1e-3
local Um = 0.3
if dim == 3 then Um = 0.45 end
local H = 0.41
local L = 0.1
local Umean2 = math.pow(2/3*Um, 2)

local ref = {}
ref.CD = 5.57953523384
ref.CL = 0.010618948146
ref.DeltaP = 0.11752016697

if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/cylinder.ugx")
	--gridName = util.GetParam("-grid", "grids/box.ugx")
	--gridName = util.GetParam("-grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_rotate_box_tri_fine.ugx")
elseif dim == 3 then
	gridName = util.GetParam("-grid", "grids/cylinder3d.ugx")
--	gridName = util.GetParam("-grid", "grids/cylinder3d_fine.ugx")
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
	local dom = Domain()
	LoadDomain(dom, gridName)
	
	
	-- NOTE: Projector creation in script-code is deprecated. Instead one should
	--		 add projectors to individual subsets directly in ProMesh.
	if     dim == 2 then 
		ProjectVerticesToSphere(dom, {0.2, 0.2}, 0.05, 0.001)
		falloffProjector = SphereProjector(MakeVec(0.2, 0.2, 0), 0.1, 0.15)
	elseif dim == 3 then 
		falloffProjector = CylinderProjector(MakeVec(0.5, 0.2, 0.0), MakeVec(0, 0, 1), 0.04, 0.1)
	end

	local projHandler = ProjectionHandler(dom:subset_handler())
	dom:set_refinement_projector(projHandler)

	projHandler:set_projector("Inner", falloffProjector)
	projHandler:set_projector("CylinderWall", falloffProjector)
	if dim == 3 then
		projHandler:set_projector("BackWall", falloffProjector)
		projHandler:set_projector("FrontWall", falloffProjector)
	end

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
	
	--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
	
	return dom
end

function CreateApproxSpace(dom, discType, vorder, porder)

	if porder == nil then porder = vorder -1 end

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder)
	
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
globalNSDisc = nil
function CreateDomainDisc(approxSpace, discType)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	globalNSDisc = NavierStokesDisc
				
	local porder = approxSpace:lfeid(dim):order()
	local vorder = approxSpace:lfeid(0):order()
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(true)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		NavierStokesDisc:set_stabilization(10)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	
	-- setup Outlet
	--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	--OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		return 4 * Um * y * (H-y) / (H*H), 0.0
	end
	function inletVel3d(x, y, z, t)
		return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0
	end
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
	--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	if dim == 2 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall")
	elseif dim == 3 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall")	
	end
	
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

function CreateSolver(approxSpace, discType)

	local base = SuperLU()
	
	local smoother = nil
	if discType == "fvcr" or discType == "fecr" then 
		smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	elseif discType == "fv1" then 
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
	end
	
	smootherDesc = util.smooth.parseParams()
--	smoother = util.smooth.create(smootherDesc)

	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	gmg:add_prolongation_post_process(AverageComponent("p"))
	transfer = StdTransfer()
	transfer:enable_p1_lagrange_optimization(false)
	gmg:set_transfer(transfer)
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))	
	end
		
	local convCheck = ConvCheck(500, 1e-11, 1e-99, true)
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, false, false))
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

	local dom = CreateDomain()
	local plots = {}	

	local function ComputeSpace(discType, p, ppress, minLev, maxLev)

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
			local domainDisc = CreateDomainDisc(approxSpace, discType)
			local solver = CreateSolver(approxSpace, discType)
		
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
				
				ComputeNonLinearSolution(u, domainDisc, solver)
				u:check_storage_type()

				local FctCmp = approxSpace:names()
				local VelCmp = {}
				
				for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
				vtkWriter = VTKOutput()
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
				meas.DeltaP.value[lev] = PEval:evaluate_global({0.15, 0.2}) - PEval:evaluate_global( {0.25, 0.2} )
		
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
	local approxSpace = CreateApproxSpace(dom, discType, vorder, porder)
	local domainDisc = CreateDomainDisc(approxSpace, discType)
	local solver = CreateSolver(approxSpace, discType)
	--solver:set_debug(GridFunctionDebugWriter(approxSpace))
			
	print(solver:config_string())
	
	local u = GridFunction(approxSpace)
	u:set(0)
	
--	ComputeNonLinearSolution(u, CreateDomainDisc(approxSpace, "fe", p), solver)
	ComputeNonLinearSolution(u, domainDisc, solver)

	local FctCmp = approxSpace:names()
	local VelCmp = {}
	
	for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
	vtkWriter = VTKOutput()
	vtkWriter:select(VelCmp, "velocity")
	vtkWriter:select("p", "pressure")
	vtkWriter:print("Cylinder", u)

	if dim == 2 then
		local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
		local C_D = 2*DL[1]/(Umean2*L)
		local C_L = 2*DL[2]/(Umean2*L)
	
		local PEval = GlobalGridFunctionNumberData(u, "p")
		local Delta_P = PEval:evaluate_global({0.15, 0.2}) - PEval:evaluate_global( {0.25, 0.2} )
	
		print("p1: "..PEval:evaluate_global({0.15, 0.2}))
		print("p2: "..PEval:evaluate_global({0.25, 0.2}))
	
		print("C_D - ref.CD: "..string.format("%.3e", C_D - ref.CD))
		print("C_L - ref.CL: "..string.format("%.3e", C_L - ref.CL))
		print("Delta_P - ref.DeltaP: "..string.format("%.3e", Delta_P - ref.DeltaP))
	end	
end
