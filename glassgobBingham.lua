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

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
order 		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
dt 			= util.GetParamNumber("-dt", 0.01)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 109)
numRefs    = util.GetParamNumber("-numRefs",    0, "number of refinements")

tAS = 0.5;
tAMR = 0.6;
tAD = 0.99;


-- Material constants
TRinne = 120.0
TUmg = 40.0
kB = 5.67e-8
--density = 2350.0
specificHeat = 1500.0
wuek = 19.0
wlf = 0.0--1280

if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/gob16x8_localrefined.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

discType = "fv1"

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    type             = " .. discType)
print("    order            = " .. order)
print("    numRefs          = "	.. numRefs)
--------------------------------------------------------------------------------
-- External Functions
--------------------------------------------------------------------------------

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

function Temperaturleitfaehigkeit(x, y, t)
	tlf = wlf/(specificHeat*density(x,y,t))
	return tlf, 0.0, 0.0, tlf
end

function density(x, y, t)
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct");
	temperature = EvalTemperature:evaluate({x,y});
	return 2375.0-0.0005*(temperature-1050.0);
	--return 2350.0
end

function VogelFulcherTammannEquation(x, y, t)
	-- calculate temperature dependent viscosity
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct");
	temperature = EvalTemperature:evaluate({x,y});
	eta = math.pow(10.0,(-2.722+4557.65/(temperature-238.541)))
	--eta = eta/density(x,y,t);
	--if eta < 0.035 or eta > 1.0 then
--		print("hier läuft was falsch: "..eta)
--	end
	return eta;
end

function contactTemperature(x, y, t)
	eps = 1e-8;
	local Eval = GlobalGridFunctionNumberData(T, "TFct")
	local T = Eval:evaluate({x,y})
	local flussRad = -0.96*kB*(math.pow(T,4)-math.pow(TUmg,4))/specificHeat
	local flussCont = -wuek * (T-TRinne)/specificHeat
	if y > Integral(1.0, u, "left")/2.0-eps or x < -Integral(1.0, u, "top")/2.0+eps or x > Integral(1.0, u, "top")/2.0-eps then		--return radiation on "top", "left" and "right"
		return flussRad
	else
		if y < -Integral(1.0, u, "left")/2.0 +eps then	--if on "bottom"
			if t < tAS+eps then return flussRad			--before Scoop
			else
				if t > tAS+eps and t < tAD+eps then return flussCont	--Scoop and 
				else return 0.0, flussRad
				end
			end
		else
			return 0.0	--return 0.0 if (x,y) is inside domain
		end
	end
end

function radiationTemperature(x, y, t)
	local Eval = GlobalGridFunctionNumberData(T, "TFct")
	local T = Eval:evaluate({x,y})
	return 0.96*kB*(math.pow(T,4)-math.pow(TUmg,4))*Integral(1.0, u, "bottom")/16.0*0.1
end

function uSol2d(x, y, t) return 0.0			end
function vSol2d(x, y, t)
	if t > tAS and t < tAMR+1e-8 then return -9.81*0.385*0.5*y*6
	else
		if t > tAMR+1e-8 and t < tAD+1e-8 then return -9.81*0.385*0.5*y
		else
			if t > 0.99+1e-8 then return -9.81*0.385*0.5*y*5
			else
				return 0.0
			end
		end
	end
end
function velSol2d(x, y, t) return uSol2d(x,y,t), vSol2d(x,y,t)  end

function sourceX(x, y, t)
	return 0.0
--	if t > 0.5-1e-8 and t < 0.6+1e-8 then
--		return -9.81*6*0.5--math.cos(30)
--	else
--		if t < 0.99+1e-8 and t > 0.6+1e-8 then
--			return -9.81*0.5--*math.cos(30)
--		else
--			if t > 0.99+1e-8 then	return -9.81*5*0.5--math.cos(30)
--			else return 0.0 end
--		end
--	end
end

function sourceY(x, y, t)
	if t > 0.5-1e-8 and t < 0.6+1e-8 then
		return -9.81*6--*0.886--math.sin(30)
	else
		if t < tAD+1e-8 and t > tAMR+1e-8 then
			return -9.81--*0.886--*math.sin(30)
		else
			if t > tAD+1e-8 then	return -9.81*6--*0.886 --*math.sin(30)
			else return 0.0 end
		end
	end
end

function source2d(x, y, t)
	return sourceX(x,y,t), sourceY(x,y,t)
end

--------------------------------------------------------------------------------
-- Loading Domain
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
	
--dom = Domain()
--LoadDomain(dom, gridName)

dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, {})

-- create Approximation Space
print("Create ApproximationSpace")
approxSpaceVel = util.ns.CreateApproxSpace(dom, discType, 1, 0)

approxSpaceTemp = ApproximationSpace(dom)	
approxSpaceTemp:add_fct("TFct", "Lagrange", 1)

u = GridFunction(approxSpaceVel)
u:set(0.0)
T = GridFunction(approxSpaceTemp)
T:set(1184.0)

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

------------------
-- Velocity
------------------

rhs = LuaUserVector("source2d")

print(approxSpaceVel:names())
NavierStokesDisc = NavierStokes(approxSpaceVel:names(), {"inner"}, discType)
NavierStokesDisc:set_source(rhs)
NavierStokesDisc:set_density("density")
NavierStokesDisc:set_laplace(false)
NavierStokesDisc:set_stokes(true)
NavierStokesDisc:set_exact_jacobian(true)	-- irrelevant
NavierStokesDisc:set_kinematic_viscosity("VogelFulcherTammannEquation");
--NavierStokesDisc:set_kinematic_viscosity(130.0/2350.0);
NavierStokesDisc:set_bingham(true)
NavierStokesDisc:set_bingham_viscosity("VogelFulcherTammannEquation");
NavierStokesDisc:set_yield_stress(0.0)
	
--upwind if available
NavierStokesDisc:set_upwind("lps")
		
-- fv1 must be stablilized
NavierStokesDisc:set_stabilization("flow", "cor")
	
-- Ausfluss soll selbst berechnet werden
-- setup Outlet
--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
--OutletDisc:add("left,right")
	
OutletDisc = DirichletBoundary()
OutletDisc:add(0.0, "p", "left,right")

InletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
InletDisc:add("top")
	
--setup Walles
--WallDisc = NavierStokesWall(NavierStokesDisc)
--WallDisc:add("bottom")

WallDisc = NavierStokesWSBCFV1(NavierStokesDisc)
WallDisc:set_sliding_factor(2.0)
WallDisc:set_sliding_limit(0.0)
WallDisc:add("bottom")

Wall1Disc = DirichletBoundary()
Wall1Disc:add(0.0, "u", "bottom")

Wall2Disc = DirichletBoundary()
Wall2Disc:add(0.0, "v", "bottom")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDiscVel = DomainDiscretization(approxSpaceVel)
domainDiscVel:add(NavierStokesDisc)
domainDiscVel:add(InletDisc)
--domainDiscVel:add(WallDisc)
--domainDiscVel:add(Wall1Disc)
domainDiscVel:add(Wall2Disc)
domainDiscVel:add(OutletDisc)

------------------
-- Temperature
------------------

--rhsTemp = LuaUserVector("contactTemperature")

elemDisc = ConvectionDiffusion("TFct", "inner", "fv1")
elemDisc:set_upwind(FullUpwind())
elemDisc:set_source("contactTemperature")
elemDisc:set_diffusion("Temperaturleitfaehigkeit")

--contactDisc = NeumannBoundary("TFct")
--contactDisc:add("contactTemperature","bottom", "inner")

--radiationDisc = NeumannBoundary("TFct")
--radiationDisc:add("radiationTemperature","top", "inner")
--radiationDisc:add("radiationTemperature","left", "inner")
--radiationDisc:add("radiationTemperature","right", "inner")

domainDiscTemp = DomainDiscretization(approxSpaceTemp)
domainDiscTemp:add(elemDisc)
--domainDiscTemp:add(contactDisc)
--domainDiscTemp:add(radiationDisc)

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

-- create time discretization
timeDiscVel = ThetaTimeStep(domainDiscVel)
timeDiscVel:set_theta(1.0) -- 1.0 is implicit euler
opVel = AssembledOperator(timeDiscVel)

opVel:init()

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

h = io.open("./csv/stokes_freeslip_densvar_viscvar_ref0_localref_noyield_nows.csv","w")
	
-- start
time = 0.0
step = 0
ly = 0.0175
lx = 0.085
flaeche = 4*lx*ly
	
FctCmp = approxSpaceVel:names()
--TempFct = approxSpaceTemp:names()

--print(approxSpaceTemp:names())
--for d = 1,#TempFct do print("TempFct["..d.."]: " .. TempFct[d]) end

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
	--Interpolate("visco", u, "VogelFulcherTammannEquation")
	--out:select(VogelFulcherTammannEquation,"visco"..dt)
	--out:select("VogelFulcherTammannEquation", "visco" .. dt)
	out:print("simpleGob", u,step,time)
	out:clear_selection()
	out:select("TFct", "temperature" .. dt)
	out:print("temperature", T, step, time)


	print(lx.."; "..ly)

	if time > 0.0-1e-8 then
		local Eval1 = GlobalGridFunctionNumberData(u, "u")

-- old --
--[[
		local summeR = 0.0
		local summeL = 0.0

		for i = 0.0, 10 do
			print(i..": "..Eval1:evaluate({lx-1e-12,-ly+1e-12+i*(ly*2.0)/10-1e-12}).."; "..Eval1:evaluate({-lx+1e-12,-ly+1e-12+i*(ly*2.0)/10-1e-12}))
			summeR = summeR + Eval1:evaluate({lx-1e-12,-ly+1e-12+i*(ly*2.0)/10-1e-12})
			summeL = summeL + Eval1:evaluate({-lx+1e-12,-ly+1e-12+i*(ly*2.0)/10-1e-12})
		end
		print("mittl. vL: " .. summeL/10)
		print("mittl. vR: " .. summeR/10)
		
		lxneu = 2*lx+dt*summeR/10-dt*summeL/10
		lyneu = flaeche/lxneu
		dlx = (lxneu)/(2*lx)
		dly = 2*ly/(lyneu)
		lx = lxneu/2
		ly = lyneu/2
		print("neue Länge: " .. 2*lx)
		print("neue Breite: " .. 2*ly)
		print("")
--]]		
	
			
-- new --
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
		
		h:write(time .. "	" ..NewLength .. "	" .. NewHeight .. "\n")
	end
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
	timeEndStep = os.clock()
	print("")
	print("Step took " .. timeEndStep-timeStartStep .. " seconds.")
end

timeEnd = os.clock()

print("Computation took " .. timeEnd-timeStart .. " seconds.")

h:close()

