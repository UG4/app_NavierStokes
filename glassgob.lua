--------------------------------------------------------------------------------
--
--   Lua - Script to compute a channel problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the channel problem
--
--   Author: Jonas Simon
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("navier_stokes_util.lua")

--------------------------------------------------------------------------------
-- Inputs
--------------------------------------------------------------------------------

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
order 		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
dt 			= util.GetParamNumber("-dt", 0.01)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 109)
numRefs    = util.GetParamNumber("-numRefs",    0, "number of refinements")

tAS = 0.5;
tAMR = 0.6;
tAD = 0.99;

-- stokes
bStokes 	= true--util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
-- bingham
bBingham 	= true--util.HasParamOption("-bingham", "If defined, bingham material is used")
--binghamViscosity 		= util.GetParamNumber("-bingvisc", 0.05)
yieldStress = util.GetParamNumber("-yieldstress", 50.0)
-- wall slip
bWallSlip 	= true--util.HasParamOption("-wallslip", "If defined, wall sliding boundary condition is used")
slidingFactor = util.GetParamNumber("-slidingFactor", -0.05)
slidingLimit = util.GetParamNumber("-slidingLimit", 0.0)
-- heat equation
TRinne = util.GetParamNumber("-TRinne", 120.0)
TUmg = util.GetParamNumber("-TUmgebung", 60.0)
kB = 5.67e-8	-- constant
specificHeat = util.GetParamNumber("-specificheat", 1500.0)
wuek = util.GetParamNumber("-heattransfer", 1280.0)
wlf = util.GetParamNumber("-conductivity", 19.0)

bFreeSlip = true

bTransVel = true

dateiname = "./csv/noVel_ys50_ws-5e-2.csv"

if 	dim == 2 then
	gridName = util.GetParam("-grid", "grids/gob16x8_localrefined_2.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

discType = "fv1"

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    type             = " .. discType)
print("    order            = " .. order)
print("    numRefs          = "	.. numRefs)
print("    bStokes          = " .. tostring(bStokes))
print("    bBingham         = " .. tostring(bBingham))
print("    yield stress     = " .. yieldStress)
print("    bWallSlip        = " .. tostring(bWallSlip))
print("    slidingLimit     = " .. slidingLimit)
print("    slidingFactor    = " .. slidingFactor)
--------------------------------------------------------------------------------
-- External Functions
--------------------------------------------------------------------------------

-- function to rescale domain
function RescaleDomain(dom, scaleX, scaleY, scaleZ)
	local sel = Selector(dom:grid())
	if scaleX == nil then print("No scaling parameter passed."); exit(); end
	if scaleY == nil then scaleY = scaleX end
	if scaleZ == nil then scaleZ = scaleY end

	
	if dom:get_dim() < 3 then scaleZ = 1 end
	if dom:get_dim() < 2 then scaleY = 1 end
	if dom:get_dim() < 1 then scaleX = 1 end

	write(">> Scaling Domain (dim: "..dom:get_dim()..") ...")	
	SelectDomainElements(sel, true, true, false, false, false)
	ScaleDomain(dom, sel, MakeVec(0, 0, 0), MakeVec(scaleX, scaleY, scaleZ)) 
	print(" done.")
end

-- calculate temperature conduction depending on density
function Temperaturleitfaehigkeit(x, y, t)
	tlf = wlf/(specificHeat*density(x,y,t))
	return tlf, 0.0, 0.0, tlf
end

-- calculate density depending on temperature
function density(x, y, t)
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct");
	temperature = EvalTemperature:evaluate({x,y});
	--return 2375.0-0.0005*(temperature-1050.0);
	return 2350.0
end

-- calculate kinematic viscosity (divided by density) depending on temperature
function VogelFulcherTammannEquationKinetic(x, y, t)
	-- calculate temperature dependent viscosity
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct");
	temperature = EvalTemperature:evaluate({x,y});
	eta = math.pow(10.0,(-2.722+4557.65/(temperature-238.541)))
	eta = eta/density(x,y,t);
	return eta;
end

-- calculate viscosity depending on temperature
function VogelFulcherTammannEquation(x, y, t)
	-- calculate temperature dependent viscosity
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct");
	temperature = EvalTemperature:evaluate({x,y});
	eta = math.pow(10.0,(-2.722+4557.65/(temperature-238.541)))
	return eta;
end

-- boundary condition for heat equation
function contactTemperature(x, y, t)
	eps = 1e-8;
	local Eval = GlobalGridFunctionNumberData(T, "TFct")
	local T = Eval:evaluate({x,y})
	local flussRad = -0.96*kB*(math.pow(T,4)-math.pow(TUmg,4))/specificHeat
	local flussCont = -wuek * (T-TRinne)/specificHeat
	if x < -Integral(1.0, u, "top")/2.0+eps or x > Integral(1.0, u, "top")/2.0-eps then		--return radiation on "top", "left" and "right"
		return flussRad
	else
		if y < -Integral(1.0, u, "left")/2.0 +eps then	--if on "bottom"
			if t < tAS+eps then return flussRad			--before Scoop
			else
				if t > tAS+eps and t < tAD+eps then return flussCont	--Scoop and 
				else return flussRad
				end
			end
		else
			if y > Integral(1.0, u, "left")/2.0-eps then	--if in "top"
				if t > tAD+eps then return flussCont
				else return flussRad
				end
			else
				return 0.0	--return 0.0 if (x,y) is inside domain
			end
		end
	end
end


function sourceX(x, y, t)
	if t < tAS+(tAMR-tAS)/2+1e-8 then
		return -9.81
	else
		if t > tAS+(tAMR-tAS)/2+1e-8 and t < tAD+1e-8 then
			return -9.81*0.5 --sin(30°) = 0.5
		else
			return -9.81
		end
	end
end


function sourceY(x, y, t)
	if t > tAS-1e-8 and t < tAMR+1e-8 then
		return -9.81*6
	else
		if t < tAD+1e-8 and t > tAMR+1e-8 then
			return -9.81*math.cos(30)
		else
			if t > tAD+1e-8 then	return 9.81*6--*0.886 --*math.sin(30)
			else return 0.0 end
		end
	end
end

function source2d(x, y, t)
	if (bTransVel) then
		return sourceX(x,y,t), sourceY(x,y,t)
	else
		return 0.0, sourceY(x,y,t)
	end
end

function initialProfile(x, y, t)
	return 1184.0-550000000.0*math.pow(y,4)
end
--------------------------------------------------------------------------------
-- Loading Domain
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
	
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, {})

-- create Approximation Space
print("Create ApproximationSpace")
approxSpaceVel = util.ns.CreateApproxSpace(dom, discType, 1, 0)

approxSpaceTemp = ApproximationSpace(dom)	
approxSpaceTemp:add_fct("TFct", "Lagrange", 1)

--------------------------------------------------------------------------------
-- Initial Conditions
--------------------------------------------------------------------------------

u = GridFunction(approxSpaceVel)
u:set(0.0)
T = GridFunction(approxSpaceTemp)
Interpolate("initialProfile", T, "TFct")
time = 0.0
step = 0

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

------------------
-- Velocity
------------------

rhs = LuaUserVector("source2d")

-- free fall
print(approxSpaceVel:names())
NavierStokesFreeFallDisc = NavierStokes(approxSpaceVel:names(), {"inner"}, discType)
NavierStokesFreeFallDisc:set_source(rhs)
NavierStokesFreeFallDisc:set_density("density")
NavierStokesFreeFallDisc:set_laplace(false)
NavierStokesFreeFallDisc:set_stokes(bStokes)
NavierStokesFreeFallDisc:set_exact_jacobian(true)	-- irrelevant
NavierStokesFreeFallDisc:set_kinematic_viscosity("VogelFulcherTammannEquationKinetic");
NavierStokesFreeFallDisc:set_bingham(bBingham)
NavierStokesFreeFallDisc:set_bingham_viscosity("VogelFulcherTammannEquation");
NavierStokesFreeFallDisc:set_yield_stress(yieldStress)
	
--upwind if available
NavierStokesFreeFallDisc:set_upwind("lps")
		
-- fv1 must be stablilized
NavierStokesFreeFallDisc:set_stabilization("flow", "cor")

-- scoop and ?
NavierStokesScoopDisc = NavierStokes(approxSpaceVel:names(), {"inner"}, discType)
NavierStokesScoopDisc:set_source(rhs)
NavierStokesScoopDisc:set_density("density")
NavierStokesScoopDisc:set_laplace(false)
NavierStokesScoopDisc:set_stokes(bStokes)
NavierStokesScoopDisc:set_exact_jacobian(true)	-- irrelevant
NavierStokesScoopDisc:set_kinematic_viscosity("VogelFulcherTammannEquationKinetic");
NavierStokesScoopDisc:set_bingham(bBingham)
NavierStokesScoopDisc:set_bingham_viscosity("VogelFulcherTammannEquation");
NavierStokesScoopDisc:set_yield_stress(yieldStress)
	
--upwind if available
NavierStokesScoopDisc:set_upwind("lps")
		
-- fv1 must be stablilized
NavierStokesScoopDisc:set_stabilization("flow", "cor")

-- deflector
NavierStokesDefDisc = NavierStokes(approxSpaceVel:names(), {"inner"}, discType)
NavierStokesDefDisc:set_source(rhs)
NavierStokesDefDisc:set_density("density")
NavierStokesDefDisc:set_laplace(false)
NavierStokesDefDisc:set_stokes(bStokes)
NavierStokesDefDisc:set_exact_jacobian(true)	-- irrelevant
NavierStokesDefDisc:set_kinematic_viscosity("VogelFulcherTammannEquationKinetic");
NavierStokesDefDisc:set_bingham(bBingham)
NavierStokesDefDisc:set_bingham_viscosity("VogelFulcherTammannEquation");
NavierStokesDefDisc:set_yield_stress(yieldStress)
	
--upwind if available
NavierStokesDefDisc:set_upwind("lps")
		
-- fv1 must be stablilized
NavierStokesDefDisc:set_stabilization("flow", "cor")

------------------
-- Temperature
------------------

-- setup boundaries
OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesFreeFallDisc)
OutletDisc:add("left,right")

BottomDisc = NavierStokesWSBCFV1(NavierStokesScoopDisc)
BottomDisc:set_sliding_factor(slidingFactor)
BottomDisc:set_sliding_limit(slidingLimit)
BottomDisc:add("bottom")

TopDefDisc = NavierStokesWSBCFV1(NavierStokesDefDisc)
TopDefDisc:set_sliding_factor(slidingFactor)
TopDefDisc:set_sliding_limit(slidingLimit)
TopDefDisc:add("top")

BottomFreeDisc = DirichletBoundary()
BottomFreeDisc:add(0.0, "p", "bottom")

TopFreeDisc = DirichletBoundary()
TopFreeDisc:add(0.0, "p", "top")

NoSlipBottomDisc = DirichletBoundary()
NoSlipBottomDisc:add(0.0, "u", "bottom")
NoSlipBottomDisc:add(0.0, "v", "bottom")

NoSlipTopDisc = DirichletBoundary()
NoSlipTopDisc:add(0.0, "u", "top")
NoSlipTopDisc:add(0.0, "v", "top")

FreeSlipBottomDisc = DirichletBoundary()
FreeSlipBottomDisc:add(0.0, "v", "bottom")

FreeSlipTopDisc = DirichletBoundary()
FreeSlipTopDisc:add(0.0, "v", "top")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDiscFreeFall = DomainDiscretization(approxSpaceVel)
domainDiscFreeFall:add(NavierStokesFreeFallDisc)
domainDiscFreeFall:add(TopFreeDisc)
domainDiscFreeFall:add(BottomFreeDisc)
domainDiscFreeFall:add(OutletDisc)

domainDiscScoop = DomainDiscretization(approxSpaceVel)
domainDiscScoop:add(NavierStokesScoopDisc)
domainDiscScoop:add(OutletDisc)
domainDiscScoop:add(TopFreeDisc)
if (bWallSlip) then
	domainDiscScoop:add(BottomDisc)
else
	if (bFreeSlip) then
		domainDiscScoop:add(FreeSlipBottomDisc)
	else
		domainDiscScoop:add(NoSlipBottomDisc)
	end
end

domainDiscDef = DomainDiscretization(approxSpaceVel)
domainDiscDef:add(NavierStokesDefDisc)
domainDiscDef:add(OutletDisc)
domainDiscDef:add(BottomFreeDisc)
if (bWallSlip) then
	domainDiscDef:add(TopDefDisc)
else
	if (bFreeSlip) then
		domainDiscDef:add(FreeSlipTopDisc)
	else
		domainDiscDef:add(NoSlipTopDisc)
	end
end


heatDisc = ConvectionDiffusion("TFct", "inner", "fv1")
heatDisc:set_upwind(FullUpwind())
heatDisc:set_source("contactTemperature")
heatDisc:set_diffusion("Temperaturleitfaehigkeit")

domainDiscTemp = DomainDiscretization(approxSpaceTemp)
domainDiscTemp:add(heatDisc)

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

-- create time discretization
timeDiscVel = ThetaTimeStep(domainDiscFreeFall)
timeDiscVel:set_theta(1.0) -- 1.0 is implicit euler
opVel = AssembledOperator(timeDiscVel)
opVel:init()

timeDiscScoop = ThetaTimeStep(domainDiscScoop)
timeDiscScoop:set_theta(1.0) -- 1.0 is implicit euler
opScoop = AssembledOperator(timeDiscScoop)
opScoop:init()

newtonSolverVel = NewtonSolver()
newtonSolverVel:set_linear_solver(LU())
newtonSolverVel:set_convergence_check(ConvCheck(500, 1e-10, 1e-99, true))

newtonSolverVel:init(opVel)

-- create time discretization
timeDiscTemp = ThetaTimeStep(domainDiscTemp)
timeDiscTemp:set_theta(1.0) -- 1.0 is implicit euler
opTemp = AssembledOperator(timeDiscTemp)

opTemp:init()

newtonSolverTemp = NewtonSolver()
newtonSolverTemp:set_linear_solver(LU())
newtonSolverTemp:set_convergence_check(ConvCheck(500, 1e-10, 1e-99, true))

newtonSolverTemp:init(opTemp)


--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

	
FctCmp = approxSpaceVel:names()

-- initialise output files
h = io.open(dateiname,"w")
out = VTKOutput()
out:clear_selection()
out:select_all(false)

-- create new grid function for old value
uOld = u:clone()
TOld = T:clone()

timeStart = os.clock()
print("Starting solver")

-- store grid functions in vector of  old solutions
solTimeSeriesVel = SolutionTimeSeries()
solTimeSeriesVel:push(uOld, time)

solTimeSeriesTemp = SolutionTimeSeries()
solTimeSeriesTemp:push(TOld, time)

for step = 1, numTimeSteps do
	timeStartStep = os.clock()
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	-- choose time step
	do_dt = dt
	
	-- setup time Disc for old solutions and timestep
	timeDiscTemp:prepare_step(solTimeSeriesTemp, do_dt)
	
	-- prepare newton solver
	if newtonSolverTemp:prepare(T) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 
	
	-- apply newton solver
	if newtonSolverTemp:apply(T) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 

	-- setup time Disc for old solutions and timestep
	timeDiscVel:prepare_step(solTimeSeriesVel, do_dt)

	-- prepare newton solver
	if newtonSolverVel:prepare(u) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 
	
	-- apply newton solver
	if newtonSolverVel:apply(u) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 

	-- update new time
	time = solTimeSeriesVel:time(0) + do_dt
		
	-- get oldest solution
	oldestSolTemp = solTimeSeriesTemp:oldest()
	oldestSolVel = solTimeSeriesVel:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSolTemp, 1.0, T)
	VecScaleAssign(oldestSolVel, 1.0, u)

	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeriesTemp:push_discard_oldest(oldestSolTemp, time)
	solTimeSeriesVel:push_discard_oldest(oldestSolVel, time)

	-- write time step data to file
	local VelCmp = {}
	for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end

	out:clear_selection()
	out:select(VelCmp, "velocity" .. dt)
	out:select("p", "pressure" .. dt)
	out:print("simpleGob", u,step,time)
	out:clear_selection()
	out:select("TFct", "temperature" .. dt)
	out:print("temperature", T, step, time)

	print("Integral Velocity (Left): " .. Integral(u, "u", "left") )
	print("Volume            (Left): " .. Integral(1.0, u, "left") )
	print("Average Velocity  (Left): " .. Integral(u, "u", "left") / Integral(1.0, u, "left") )

	print("Integral Velocity (Right): " .. Integral(u, "u", "right") )
	print("Volume            (Right): " .. Integral(1.0, u, "right") )
	print("Average Velocity  (Right): " .. Integral(u, "u", "right") / Integral(1.0, u, "right") )

	local AvgVelLeft = Integral(u, "u", "left") / Integral(1.0, u, "left");
	local AvgVelRight = Integral(u, "u", "right") / Integral(1.0, u, "right");

	local Height = Integral(1.0, u, "left")
	local Length = Integral(1.0, u, "bottom")
		
	local NewLength = Length + dt*(AvgVelRight - AvgVelLeft)
	local NewHeight = (Length*Height) / NewLength
		
	local dlx = NewLength / Length
	local dly = NewHeight / Height 

	print("Scale: x: " .. dlx .. ", y: " .. dly)
	print("Neue Länge: " .. NewLength)
	print("Neue Höhe: " .. NewHeight)
	print(" ")

	RescaleDomain(dom, dlx, dly, 0.0)
		
	h:write(time .. "	" ..NewLength .. "	" .. AvgVelLeft .. "	" .. AvgVelRight .. "\n")
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
	timeEndStep = os.clock()
	print("")
	print("Step took " .. timeEndStep-timeStartStep .. " seconds.")
	print("")

	-- -0.01 because we change at end of step i-1 for step i
	if time > tAS-0.01-1e-8 and time < tAS-0.01+1e-8 then
		timeDiscVel = ThetaTimeStep(domainDiscScoop)
		timeDiscVel:set_theta(1.0) -- 1.0 is implicit euler
		opVel = AssembledOperator(timeDiscVel)
		opVel:init()

		newtonSolverVel:init(opVel)
	end
	if time > tAD-1e-8 and time < tAD+1e-8 then 
		timeDiscVel = ThetaTimeStep(domainDiscDef)
		timeDiscVel:set_theta(1.0)
		opVel = AssembledOperator(timeDiscVel)
		opVel:init()

		newtonSolverVel:init(opVel)
	end
end

timeEnd = os.clock()

print("Computation took " .. timeEnd-timeStart .. " seconds.")

h:close()

