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

-- Material constants
TRinne = 120.0
kB = 5.67e-8
density = 2350.0
specificHeat = 1500.0

if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/gob16x8.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

discType = "fv1"

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    type             = " .. discType)
print("    order            = " .. order)

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

function VogelFulcherTammannEquation(x, y, t)
	-- calculate temperature dependent viscosity
	EvalTemperature = GlobalGridFunctionNumberData(T, "TFct")
	temperature = EvalTemperature:evaluate({x,y})
	eta = math.pow(10.0,(-2.722+4557.65/(temperature-238.541)))
	return eta/density
end

function temperatureBC(x, y, t)
	local Eval = GlobalGridFunctionNumberData(T, "TFct")
	local T = Eval:evaluate({x,y})
	local fluss = -0.96*kB*(math.pow(T,4)-math.pow(TRinne,4))*(0.17/16.0)
	if y < -0.0175+1e-5 then
		return fluss*0.1
	else
		if x < -(0.17/2)+1e-5 or x > (0.17/2)-1e-5 or y > 0.0175-1e-5 then
			return fluss*0.5
		else
			return 0.0
		end
	end
end

function uSol2d(x, y, t) return 0.0			end
function vSol2d(x, y, t)
	if t > 0.5 and t < 0.6+1e-8 then return -9.81*0.385*0.5*y*6
	else
		if t > 0.6+1e-8 and t < 0.99+1e-8 then return -9.81*0.385*0.5*y
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
		return -9.81*6*0.886*density--math.sin(30)
	else
		if t < 0.99+1e-8 and t > 0.6+1e-8 then
			return -9.81*0.886*density--*math.sin(30)
		else
			if t > 0.99+1e-8 then	return -9.81*5*0.886*density--*math.sin(30)
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
	
dom = Domain()
LoadDomain(dom, gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpaceVel = util.ns.CreateApproxSpace(dom, discType, 1, 0)

approxSpaceTemp = ApproximationSpace(dom)	
approxSpaceTemp:add_fct("TFct", "Lagrange", 1)

u = GridFunction(approxSpaceVel)
u:set(0)
T = GridFunction(approxSpaceTemp)
T:set(1184.0)
--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

rhs = LuaUserVector("source2d")
print(rhs)

print(approxSpaceVel:names())
NavierStokesDisc = NavierStokes(approxSpaceVel:names(), {"inner"}, discType)
NavierStokesDisc:set_source(rhs)
NavierStokesDisc:set_density(2350.0)
NavierStokesDisc:set_stokes(true)
NavierStokesDisc:set_exact_jacobian(true)	--??
NavierStokesDisc:set_kinematic_viscosity("VogelFulcherTammannEquation");
	
--upwind if available
NavierStokesDisc:set_upwind("full")
		
-- fv1 must be stablilized
NavierStokesDisc:set_stabilization("flow", "cor")
	
-- Ausfluss soll selbst berechnet werden
-- setup Outlet
OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDisc:add("left,right")
	

--OutletDisc = DirichletBoundary()
--OutletDisc:add(0.0, "p", "left,right")

-- setup Inlet
InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("velSol"..dim.."d", "top")
	
--setup Walles
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("bottom")
Wall1Disc = DirichletBoundary()
Wall1Disc:add(0.0, "u", "bottom")

Wall2Disc = DirichletBoundary()
Wall2Disc:add(0.0, "v", "bottom")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDiscVel = DomainDiscretization(approxSpaceVel)
domainDiscVel:add(NavierStokesDisc)
--domainDiscVel:add(InletDisc)
--domainDiscVel:add(WallDisc)
domainDiscVel:add(Wall1Disc)
domainDiscVel:add(Wall2Disc)
domainDiscVel:add(OutletDisc)
	

elemDisc = ConvectionDiffusion("TFct", "inner", "fv1")
elemDisc:set_upwind(FullUpwind())
elemDisc:set_source("temperatureBC")
elemDisc:set_diffusion(0.00028)

domainDiscTemp = DomainDiscretization(approxSpaceTemp)
domainDiscTemp:add(elemDisc)

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
newtonSolverVel:set_convergence_check(ConvCheck(500, 1e-11, 1e-99, true))

newtonSolverVel:init(opVel)

-- create time discretization
timeDiscTemp = ThetaTimeStep(domainDiscTemp)
timeDiscTemp:set_theta(1.0) -- 1.0 is implicit euler
opTemp = AssembledOperator(timeDiscTemp)

opTemp:init()

newtonSolverTemp = NewtonSolver()
newtonSolverTemp:set_linear_solver(LU())
newtonSolverTemp:set_convergence_check(ConvCheck(500, 1e-11, 1e-99, true))

newtonSolverTemp:init(opTemp)


--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

h = io.open("laengen_freeslip.csv","w")
	
-- start
time = 0.0
step = 0
ly = 0.0175
lx = 0.085
flaeche = 4*lx*ly
	
FctCmp = approxSpaceVel:names()

out = VTKOutput()
out:clear_selection()
out:select_all(false)

-- create new grid function for old value
uOld = u:clone()
TOld = T:clone()

timeStart = os.clock()
print("Starting solver")

-- store grid function in vector of  old solutions
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
	if newtonSolverTemp:prepare(u) == false then 
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

	print(lx.."; "..ly)

	if time > 0.0-1e-8 then
		local Eval1 = GlobalGridFunctionNumberData(u, "u")

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
		RescaleDomain(dom, dlx, dly, 0.0)
	end
	h:write(time .. "	" .. 2*lx .. "	" .. 2*ly .. "\n")
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
	timeEndStep = os.clock()
	print("")
	print("Step took " .. timeEndStep-timeStartStep .. " seconds.")
end

timeEnd = os.clock()

print("Computation took " .. timeEnd-timeStart .. " seconds.")

h:close()

