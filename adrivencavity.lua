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
exactJacFactor = util.GetParamNumber("-exactjac", 0)
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "no", "Upwind type")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "COR", "Diffusion length type")
bPSep       = util.HasParamOption("-psep", "If defined, pressure separation used")
bPLin       = util.HasParamOption("-linp", "If defined, pressure gradient is used")
bPLinDefect   = util.HasParamOption("-linpdefect", "If defined, pressure gradient is used only in defect")
bNoUpwindInDefect = util.HasParamOption("-noupdefect", "If defined, no upwind is used in defect")
bLinUpwindInDefect = util.HasParamOption("-linupdefect", "If defined, linear upwind is used in defect")
linred      = util.GetParam("-linred", 1e-1 , "Linear reduction")
nlintol     = util.GetParam("-nlintol", 1e-5, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")
bNoLineSearch  = util.HasParamOption("-noline", "If defined, no line search is used")
boolStat = util.GetParamNumber("-stat", 1)
elemType = util.GetParam("-elem", "quad")
dt = util.GetParamNumber("-dt", 0.1)
reynoldsNr = util.GetParamNumber("-R",100)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)
numAdaptRefs  = util.GetParamNumber("-numAdaptRefs",0)

if type==fv1 then
	InitUG(dim, AlgebraType("CPU", dim+1));
else
	InitUG(dim, AlgebraType("CPU", 1));
end

-- undo fvcr only options if type is not fvcr
if upwind == "linear" then
	if type~="fvcr" then
		print("Upwind type '"..upwind.."' only supported for fvcr discretization.")
		upwind = "no"
	else
		upwind = "full"
		bLinearUpwind = true
		bLinUpwindInDefect = true
	end
else
	bLinearUpwind = false
end

if bPSep == true then
	if type~="fvcr" then
		print("Warning: pressure separation only supported for fvcr discretization.")
		bPSep=false
	end
end

if bPLin == true then
	if type~="fvcr" then
		print("Warning: pressure gradient defect modification only supported for fvcr discretization.")
		bPLin=false
	end
end


if upwind == "linear" then
	bLinUpwindInDefect=true
end
if bPLin == true then
	bPLinDefect=true
end
if bPLinDefect==false then
	bPLin=false
end

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "grids/unit_square_01_tri_unstruct_4bnd.ugx")
	else 
		gridName = util.GetParam("-grid", "grids/dc_quads.ugx")
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2_four_bnd.ugx")	
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

if bNoUpwindInDefect == true then
	if type~="fvcr" then
		print("Warning: no upwind in defect modification only supported for fvcr discretization.")
		bNoUpwindInDefect=false
	else
		upwind = "full"
	end
end
if bLinUpwindInDefect == true then
	if type~="fvcr" then
		print("Warning: linear upwind in defect modification only supported for fvcr discretization.")
		bLinUpwindInDefect=false
	else
		upwind = "full"
	end
end


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
print("    exact jacob. factor = " .. exactJacFactor)
print("    peclet blend        = " .. tostring(bPecletBlend))
print("    upwind              = " .. upwind)
if type=="fv1" then
	print("    pac upwind          = " .. tostring(bPac))
	print("    stab                = " .. stab)
	print("    diffLength          = " .. diffLength)
end
if type=="fvcr" then
	print("    pressure separation = " .. tostring(bPSep))
	print("    no upwind in defect = " .. tostring(bNoUpwindInDefect))
	if bLinUpwind==true then
		print("    linear upwind         = " .. tostring(bLinUpwind))
	else
		print("    linear upwind in def  = " .. tostring(bLinUpwindInDefect))
	end
	if bPLin==true then
		print("    linear pressure       = " .. tostring(bPLin))
	else
		print("    lin pressure in def   = " .. tostring(bPLinDefect))
	end
end
print("    no line search      = " .. tostring(bNoLineSearch))
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

if numAdaptRefs>0 then
	refiner = HangingNodeDomainRefiner(dom);
end

IdentifySubsets(dom,"Left","Right")

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

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", type)
NavierStokesDisc:set_exact_jacobian(exactJacFactor)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(1.0/reynoldsNr)




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



--------------------------------------------------------------------------------
-- Boundary conditions and constraints
--------------------------------------------------------------------------------

InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add({1,0}, "Top")

WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("Bottom")

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
	timeDisc:set_theta(0.5) -- 1.0 is implicit euler
	op = AssembledOperator(timeDisc)
end
op:init()

u = GridFunction(approxSpace)

if bNoUpwindInDefect == true then
	NavierStokesDisc:set_defect_upwind(false)
end
if type=="fvcr" then
	tBefore = os.clock()
--	OrderCRCuthillMcKee(approxSpace,u,true)
--	CROrderCuthillMcKee(approxSpace,u,true,false,false,true)
	-- CROrderSloan(approxSpace,u,false,false,true)
	-- CROrderKing(approxSpace,u,true,false,false,true)
--	CROrderMinimumDegree(approxSpace,u,true,false,true)
	-- OrderLex(approxSpace, "lr")
	tAfter = os.clock()
	tOrder = tAfter-tBefore
	print("Ordering took " .. tAfter-tBefore .. " seconds.")
else
	if type=="fv1" then
		OrderCuthillMcKee(approxSpace,true)
	end
end

function StartValue_u2d(x,y,t) return math.sin(2*math.pi*x)*y*y*y*math.exp(x) end
function StartValue_v2d(x,y,t) return 0 end

Interpolate("StartValue_u"..dim.."d", u, "u")
Interpolate("StartValue_v"..dim.."d", u, "v")

if (bLinearUpwind==true)or(bPLin==true)or(bLinUpwindInDefect==true) then
	domainDisc:add(DiscConstraintFVCR(u,bLinUpwindInDefect,bLinearUpwind,bPLinDefect,bPLin,false))
end

egsSolver = LinearSolver()
egsSolver:set_preconditioner(ElementGaussSeidel("vertex"))
-- egsSolver:set_preconditioner(Vanka())
-- egsSolver:set_preconditioner(ILUT(1e-16))
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
	ilutsmoother:set_threshold(1e-4)
--	ilutsmoother:set_info(true)
	smoother=ilutsmoother
elseif type=="fvcr" then 
	smoother=Vanka()
	smoother=CRILUT(1e-1,1e-3,true)
elseif type=="fv" or type=="fe" then 
	smoother=ElementGaussSeidel()
end
smoother:set_damp(0.8)

if type=="fv1" then 
	basePre = ILUT()
	basePre:set_threshold(1e-7)
elseif type=="fvcr" then 
	basePre = CRILUT(1e-16,false)
elseif type=="fv" or type=="fe" then 
	basePre=ElementGaussSeidel()
end
baseSolver = LinearSolver()
baseSolver:set_preconditioner(basePre)
baseSolver:set_convergence_check(ConvCheck(10000, lintol*0.1,linred*0.1,true))

if type=="fv1" then   
--  for ilu 2 is sufficient
	numSmooth=2
else
--  for Vanka type smoothers there should be more smoothing steps on higher levels
	numSmooth=2*numRefs
	numSmooth=10
end

VankaSolver = LinearSolver()
VankaSolver:set_preconditioner(Vanka())
VankaSolver:set_convergence_check(ConvCheck(100, lintol, 1, false))

ilutSolver  = LinearSolver()
ilutSolver:set_preconditioner(CRILUT(1e-3))
ilutSolver:set_convergence_check(ConvCheck(10, lintol, 0.9999999, true))
-- dbgWriter = GridFunctionDebugWriter(approxSpace)
-- ilutSolver:set_debug(dbgWriter)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(2)
gmg:set_base_solver(LU())
gmg:set_base_solver(VankaSolver)
-- gmg:set_base_solver(ilutSolver)
gmg:set_smoother(smoother)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(numSmooth)
gmg:set_num_postsmooth(numSmooth)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_emulate_full_refined_grid(true)
-- gmg:set_smooth_on_surface_rim(true)
gmg:set_rap(false)
-- gmg:set_debug(dbgWriter)
gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(100, lintol, linred, true))
-- gmgSolver:set_debug(dbgWriter)

BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(smoother)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, lintol, linred, true))

productPre = LinearIteratorProduct()
productPre:add_iterator(gmg)
productPre:add_iterator(CRILUT(1e-1,1e-3),4)
productSolver = LinearSolver()
productSolver:set_preconditioner(productPre)
productSolver:set_convergence_check(ConvCheck(100000, lintol, linred, true))

if type=="fv1" then 
	solver=gmgSolver
elseif type=="fvcr" then 
	solver=ilutSolver
	solver=gmgSolver
	solver=productSolver
--	solver=BiCGStabSolver
--	solver=egsSolver
elseif type=="fv" or type=="fe" then 
	solver=productSolver
--	solver=ilutSolver
--	solver=LU()
	solver=gmgSolver
--	solver=BiCGStabSolver
end

-- solver = VankaSolver

newtonSolver = NewtonSolver(op)
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(ConvCheck(250, nlintol, nlinred, true))
if bNoLineSearch==false then
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true))
end

-- newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- SaveVectorForConnectionViewer(u, "StartSolution.vec")

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

targetVertexNumber = 10000
alpha = 0.5
coarsenTol = 0.00000000001
refTol = 1e-2

-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()
if boolStat==1 then

	newtonSolver:init(op)
		if newtonSolver:prepare(u) == false then
		print ("Newton solver prepare failed.") exit()
	end

--	SaveVectorForConnectionViewer(u, "StartSolution.vec")

	if newtonSolver:apply(u) == false then
		print ("Newton solver apply failed.") exit()
	end
	
	-- pressure separation (only fvcr)
	if bPSep==true then
		pGradSource=SeparatedPressureSource(approxSpace,u)
		NavierStokesDisc:set_source(pGradSource)
		pGradSource:update()
		
		dcevaluation(u,reynoldsNr)
		
		function zero(x,y,t) return 0 end
		Interpolate("zero", u, "p")
			
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step..".") exit() 
		end 
		
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			print ("Newton solver failed at step "..step..".") exit() 
		end 
	end
	
	out:print("DrivenCavity", u)
else

	-- store grid function in vector of  old solutions
	solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(uOld, time)
	
	out:print("TimeDrivenCavity", u,0,0)

	for step = 1, numTimeSteps do
		print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

		-- choose time step
		do_dt = dt
	
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
	
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step..".") exit() 
		end 
	
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			print ("Newton solver failed at step "..step..".") exit() 
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

		print("++++++ TIMESTEP " .. step .. "  END ++++++")
		
		-- plot solution
		out:print("TimeDrivenCavity", u,step,time)
		
		maxLvl = numRefs+numAdaptRefs
		
		
		
		for adaptStep = 1, numAdaptRefs do
			oldVertexNumber = ParallelSum(dom:grid():num_vertices())
			
			-- 6. estimate error and mark
			if (adaptStep==1) then
				print("vertex number " .. ParallelSum(dom:grid():num_vertices()))
				print("coarsen tolerance = " .. coarsenTol) 
				print("refinement tolerance = " .. refTol) 
			end
		
			--	mark for coarsening
			if(coarsenTol > 0) then
--				print("Estimating error")
--				refiner:clear_marks()
--				MarkForAdaption_AbsoluteGradientIndicator(refiner, u, "p", -1, coarsenTol, 0, maxLvl)
--				print("Coarsening")
--				refiner:coarsen()
			end
			
			vertexNumber = ParallelSum(dom:grid():num_vertices())
			
			print("vertex number after coarsening " .. vertexNumber)
			
			--	mark for refinement
			--refiner:set_adjusted_marks_debug_filename("adjmarks-"..step..".ugx")
			if(refTol >= 0) then
				print("Estimating error")
				refiner:clear_marks()
				MarkForAdaption_AbsoluteGradientIndicator(refiner, u, "p", refTol, -1, numRefs, maxLvl)
--				MarkForAdaption_GradientIndicator(refiner, u, "p", refTol, 0, maxLvl)
--				MarkForAdaption_GradientJumpIndicator(refiner, u, "c", refTol, desc.refFrac, 0, maxLvl)
				print("Elements marked.")
			
				-- 7. refine
				print("refining")
				local nodesBefore = ParallelSum(dom:grid():num_vertices())
				refiner:refine() 
				local nodesAfter = ParallelSum(dom:grid():num_vertices())
			
				if(nodesAfter == nodesBefore) then
					print("No new elements have been introduced during refinement.")
					print("No more refinements will be attempted!")
					breakAtEndOfLoop = true
				else
					print("Grid refined.")
				end	
			end
			
			MarkForRefinement_VerticesInSphere(dom, refiner, MakeVec(-0.5, 0) , 0.6 )
			refiner:refine()
			--SaveGridHierarchyTransformed(dom:grid(), "gh-"..step..".ugx", 0.1)
			-- 8. defragment approximation space
			
			vertexNumber = ParallelSum(dom:grid():num_vertices())
			
			print("vertex number after refinement " .. vertexNumber)
			
			ratio=oldVertexNumber/vertexNumber
			neededRatio = vertexNumber/targetVertexNumber
			beta = 0
			if (ratio<1) and (neededRatio>1) then
				beta = alpha
			end
			if (ratio>1) and (neededRatio<1) then
				beta = -alpha
			end
			coarsenTol = coarsenTol*(1+beta)
			refTol = refTol*(1+beta)
			print("new coarsen tolerance = " .. coarsenTol) 
			print("new refinement tolerance = " .. refTol) 
			
		end	
		if numAdaptRefs>0 then
			approxSpace:print_statistic()
			print("Statistic done")
		end
	end
end
tAfter = os.clock()
dcevaluation(u,reynoldsNr)
newtonSolver:print_average_convergence()
print("Computation took " .. tAfter-tBefore .. " seconds.")
print("Computation + ordering took " .. tAfter-tBefore+tOrder .. " seconds.")
