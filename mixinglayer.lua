-------------------------------------------------------------------------------------------------------
--
--  Mixing layer problem (test case for turbulence)
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "fvcr" or "stabil"
discType = util.GetParam("-type", "fvcr")
elemType = util.GetParam("-elem", "quads")

InitUG(dim, AlgebraType("CPU", 1))

-- Setup from John LES book (p. 200) computation on [-1 1] square
sigma0=1.0/14.0
winf=1
cnoise=0.001
viscosity=1/140000
timeUnit=sigma0/winf

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "grids/unit_square_01_tri_4bnd.ugx")
	else
		gridName = util.GetParam("-grid", "grids/unit_square_01_quads_2x2_4bnd.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting.") exit() end

undefined    = -3458789.116

dt = util.GetParamNumber("-dt", undefined)
dtTimeUnit = util.GetParamNumber("-dtScale", 0.1)*timeUnit
timeMethod = util.GetParam("-timeMethod","cn")
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 100)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",3)
turbViscMethod = util.GetParam("-turbulenceModel","dyn")
modellconstant = util.GetParamNumber("-c",0.1)
bLaplace 	= util.HasParamOption("-laplace", "If defined, only laplace term used")
exJacFactor = util.GetParamNumber("-exactjac", 0)
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type")
nolimit      = util.HasParamOption("-nolimit", "If defined, no limiter is used in linear upwind")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "COR", "Diffusion length type")
bPLin       = util.HasParamOption("-linp", "If defined, pressure gradient is used")
bPLinDefect   = util.HasParamOption("-linpdefect", "If defined, pressure gradient is used only in defect")
bNoUpwindInDefect = util.HasParamOption("-noupdefect", "If defined, no upwind is used in defect")
bLinUpwindInDefect = util.HasParamOption("-linupdefect", "If defined, linear upwind is used in defect")
graddivFactor = util.GetParamNumber("-graddiv", 0)
linred      = util.GetParam("-linred", 1e-1 , "Linear reduction")
nlintol     = util.GetParam("-nlintol", 1e-6, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")
bNoLineSearch  = util.HasParamOption("-noline", "If defined, no line search is used")
bSave          = util.HasParamOption("-save", "If defined solution vector is safed after every step")
outputFactor     = util.GetParam("-output", 1, "output every ... steps")
tsOffset       = util.GetParamNumber("-tsOffset",0) 
startTime       = util.GetParamNumber("-starttime",undefined) 

if dt == undefined then
	dt = dtTimeUnit
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
if (bNoUpwindInDefect == true) or (bLinUpwindInDefect == true) then
	upwind = "full"
end

if nolimit==true then
	bLimit=false
else
	bLimit=true
end

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    dt           = " .. dt)
print("    numTimeSteps = " .. numTimeSteps)
print("    time stepping method = " .. timeMethod)
print("    turbulence model    = " .. turbViscMethod)
print("    grid       	       = " .. gridName)
print("    laplace             = " .. tostring(bLaplace))
print("    grad-div factor     = " .. graddivFactor)
print("    exact jacob. factor = " .. exJacFactor)
print("    peclet blend        = " .. tostring(bPecletBlend))
print("    upwind              = " .. upwind)
print("    no upwind in defect = " .. tostring(bNoUpwindInDefect))
if bLinUpwind==true then
	print("    linear upwind         = " .. tostring(bLinUpwind))
else
	print("    linear upwind in def  = " .. tostring(bLinUpwindInDefect))
end
if bLinUpwind==true or bLinUpwindInDefect==true then
	print("    limiter               = " .. tostring(bLimit))
end
if bPLin==true then
	print("    linear pressure       = " .. tostring(bPLin))
else
	print("    lin pressure in def   = " .. tostring(bPLinDefect))
end
print("    no line search      = " .. tostring(bNoLineSearch))
print("    linear reduction    = " .. linred)
print("    linear tolerance    = " .. lintol)
print("    nonlinear reduction = " .. nlinred)
print("    nonlinear tolerance = " .. nlintol)

if upwind == "linear" then
	upwind = "full"
	bLinearUpwind = true
else
	bLinearUpwind = false
end

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
IdentifySubsets(dom,"Left","Right")

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
approxSpace:add_fct("p", "piecewise-constant") 

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

approxSpaceVorticity = ApproximationSpace(dom)
approxSpaceVorticity:add_fct("c", "Crouzeix-Raviart", 1)

vort = GridFunction(approxSpaceVorticity)

u = GridFunction(approxSpace)
OrderCRCuthillMcKee(approxSpace,u,true)
-- OrderLex(approxSpace, "lr")

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", "fvcr")

-- set parameters
NavierStokesDisc:set_upwind(upwind)
NavierStokesDisc:set_peclet_blend(bPecletBlend)
NavierStokesDisc:set_exact_jacobian(exJacFactor)
NavierStokesDisc:set_grad_div(graddivFactor)
NavierStokesDisc:set_stokes(false)
NavierStokesDisc:set_laplace(bLaplace)

----------------------------------
----------------------------------
-- Viscosity Data
----------------------------------
----------------------------------

if turbViscMethod=="no" then
	NavierStokesDisc:set_kinematic_viscosity(viscosity)
else
	if turbViscMethod=="dyn" then
		viscosityData = CRDynamicTurbViscData(approxSpace,u)
	end
	if turbViscMethod=="sma" then
		viscosityData = CRSmagorinskyTurbViscData(approxSpace,u,modellconstant)
	end	
	viscosityData:set_turbulence_zero_bnd("Top,Bottom")
	viscosityData:set_kinematic_viscosity(viscosity)
	NavierStokesDisc:set_kinematic_viscosity(viscosityData)
end


------------------------------------------
------------------------------------------
-- Boundary conditions
------------------------------------------
------------------------------------------

-- OutletDiscTop = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDiscTop = CRNavierStokesSymBC(NavierStokesDisc)
OutletDiscTop:add("Top")
-- OutletDiscBottom = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDiscBottom = CRNavierStokesSymBC(NavierStokesDisc)
OutletDiscBottom:add("Bottom")

------------------------------------------
------------------------------------------
-- Set up discretization and constraints
------------------------------------------
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(OutletDiscTop)
domainDisc:add(OutletDiscBottom)

if (bLinearUpwind==true)or(bLinUpwindInDefect==true)or(bPLin==true) then
	-- last three parameters: adaptivity boolean, gradient limiter boolean, boundaries where full upwind/constant pressure is used 
	domainDisc:add(DiscConstraintFVCR(u,bLinUpwindInDefect,bLinearUpwind,bPLinDefect,bPLin,false,bLimit,"Top,Bottom"))
end

-- create operator from discretization

-- create time discretization
if timeMethod=="cn" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.5) -- Crank-Nicolson method
end
if timeMethod=="euler" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1) -- implicit Euler
end	
if timeMethod=="fracstep" then
	timeDisc = ThetaTimeStep(domainDisc,"FracStep")
end
if timeMethod=="alex" then
	timeDisc = ThetaTimeStep(domainDisc, "Alexander")
end

op = AssembledOperator(timeDisc)

op:init()

function StartValue_u2d(x,y,t) 
	return winf*math.tanh(2*y/sigma0)+cnoise*winf*(-8*y*math.exp(-(2*y/sigma0)*(2*y/sigma0))*(math.cos(8*math.pi*x)+math.cos(20*math.pi*x)))/sigma0/sigma0
end

function StartValue_v2d(x,y,t) 
	return cnoise*winf*(-math.exp(-4*y*y/sigma0/sigma0)*(-8*math.sin(8*math.pi*x)*math.pi-20*math.sin(20*math.pi*x)*math.pi))
end

-- filter value for start solution
h = 1
for i=1,numRefs do
	h = 0.5 * h
end

function StartValueFiltered_u2d(x,y,t) 
	return 0.25*( StartValue_u2d(x+h,y)+StartValue_u2d(x-h,y)+StartValue_u2d(x,y+h)+StartValue_u2d(x,y-h) )
end

function StartValueFiltered_v2d(x,y,t) 
	return 0.25*( StartValue_v2d(x+h,y)+StartValue_v2d(x-h,y)+StartValue_v2d(x,y+h)+StartValue_v2d(x,y-h) )
end

function StartValue_p2d(x,y,t) return 0 end

-- if offset is > 0 use current solution as start solution
if tsOffset==0 then
	Interpolate("StartValueFiltered_u2d", u, "u")
	Interpolate("StartValueFiltered_v2d", u, "v")
	Interpolate("StartValue_p2d", u, "p")
else
	LoadVector(u,"currentSolution.vec")
end

vanka = Vanka()
vanka:set_damp(0.9)

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, false))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-5)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(false)

CRILUT = CRILUT()
CRILUT:set_threshold(1e-0,1e-1)
CRILUT:set_damp(0.95)
-- CRILUT:set_info(true)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(CRILUT)
ilutSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, false))

vankaBase = LinearSolver()
vankaBase:set_preconditioner(Vanka())
vankaBase:set_convergence_check(baseConvCheck)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(ilutSolver)
gmg:set_smoother(CRILUT)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
-- gmg:set_damp(MinimalResiduumDamping())
-- gmg:set_damp(0.8)
-- gmg:set_damp(MinimalEnergyDamping())

-- gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(vanka)
-- BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-5, 1e-1, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, lintol, linred, true))

-- choose a solver
solver = BiCGStabSolver
solver = vankaSolver
solver = gmgSolver
-- solver = ilutSolver

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10000)
newtonConvCheck:set_minimum_defect(nlintol)
newtonConvCheck:set_reduction(nlinred)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.7)
newtonLineSearch:set_accept_best(false)

dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
if bNoLineSearch==false then
	newtonSolver:set_line_search(newtonLineSearch)
end
-- newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)
if turbViscMethod~="no" then
	newtonSolver:add_step_update(viscosityData)
end

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end

-- if newtonSolver:apply(u) == false then
--	 print ("Newton solver apply failed.") exit()
-- end

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- start
step = 0
if startTime==undefined then
	time = tsOffset * dt
else
	time = startTime
end

if tsOffset==0 then
	-- compute initial vorticity
	vorticity(vort,u)

	-- write start solution
	print("Writing start values")
	out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	out:select_element("u,v", "velocity")
	out:select_element("u", "u")
	out:select_element("v", "v")
	out:select_element("p", "p")
	out:print("MixingLayer", u,0,0)

	outv = VTKOutput()
	outv:select_element("c","c")
	outv:print("vorticity", vort,0,0)
end

-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

function zero(x,y,t) return 0 end

--clearFile("kineticEnergy.m")
-- compute kinetic energy
ke=kineticEnergy(u)
--writeNumbers("kineticEnergy.m",1,0,ke)

for step = 1 + tsOffset, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	-- choose time step
	do_dt = dt
	
	for stage = 1, timeDisc:num_stages() do
	
		if timeDisc:num_stages() > 1 then
			timeDisc:set_stage(stage)
			print("      +++ STAGE " .. stage .. " BEGIN ++++++")
		end
	
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
	
	end
	
	print ("Time units = ".. time/timeUnit)
	
	-- compute CFL number 
	cflNumber(u,do_dt)
	
	-- compute kinetic energy
	ke=kineticEnergy(u)
	--writeNumbers("kineticEnergy.m",step+1,time,ke)
	
	if step % outputFactor == 0 then
	
		-- compute vorticity
		vort:set(0)
		vorticity(vort,u)
	
		if (bSave) then
			SaveVectorForConnectionViewer(u,"currentSolution.vec")
		end
	
		out = VTKOutput()
		out:clear_selection()
		out:select_all(false)
		out:select_element("u,v", "velocity")
		out:select_element("u", "u")
		out:select_element("v", "v")
		out:select_element("p", "p")
		out:print("MixingLayer", u,step,time)
		outv = VTKOutput()
		outv:select_element("c","c")
		outv:print("vorticity", vort,step,time)
		print(" ")
		print("++++++ TIMESTEP " .. step .. "  END ++++++")
		
	end
end

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.")

-- plot solution

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("DCSolution", u)

print("done.")
