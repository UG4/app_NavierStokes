--------------------------------------------------------------------------------
--
--   Lua - Script to test AddFunctionValuesToGridCoordinatesP1
--
--	This script sets defines a constant source and resizes the mesh pointwise
--  after every timestep. No other calculation is performed.
--
--   Author: Jonas Simon
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("navier_stokes_util.lua")

-- parameter
dim 	  		= util.GetParamNumber("-dim", 2, "world dimension")
order 	  		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
dt 			 	= util.GetParamNumber("-dt", 0.5)
numTimeSteps 	= util.GetParamNumber("-numTimeSteps", 10)

discType = "fv1"

gridName = util.GetParam("-grid", "grids/gob16x8_localrefined_2.ugx")


-- functions
function RescaleDomainPointwise(dom)
	AddFunctionValuesToGridCoordinatesP1(u, "u", 0, dt)
	AddFunctionValuesToGridCoordinatesP1(u, "v", 1, dt)
end

function xVelocity(x,y,t)
	return x/Integral(1.0, u, "bottom")*0.1
end

function yVelocity(x,y,t)
	return y/Integral(1.0, u, "left")*0.1
end

function source2d(x,y,t)
	return xVelocity(x,y,t), yVelocity(x,y,t)
end

-- initialisation
InitUG(dim, AlgebraType("CPU", 1))

dom = util.CreateAndDistributeDomain(gridName, 0, 0, {})

-- create Approximation Space
approxSpace = util.ns.CreateApproxSpace(dom, discType, 1, 0)

-- 
u = GridFunction(approxSpace)
--u:set(0.0)
Interpolate(xVelocity, u, "u")
Interpolate(yVelocity, u, "v")

rhs = LuaUserVector("source2d")

Disc = NavierStokes(approxSpace:names(), {"inner"}, discType)
Disc:set_source(rhs)
Disc:set_stokes(true)

domainDisc = DomainDiscretization(approxSpace)

timeDisc = ThetaTimeStep(domainDisc)
--timeDisc:set_theta(1.0) -- 1.0 is implicit euler
op = AssembledOperator(timeDisc)
op:init()

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(LU())
newtonSolver:set_convergence_check(ConvCheck(500, 1e-10, 1e-99, true))

newtonSolver:init(op)

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------
time = 0.0

FctCmp = approxSpace:names()

local VelCmp = {}
for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select(VelCmp, "velocity at begin" .. dt)
out:print("test", u,0,0.0)

-- create new grid function for old value
uOld = u:clone()

-- store grid functions in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, numTimeSteps do
	-- choose time step
	do_dt = dt

	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)

	-- prepare newton solver
	if newtonSolver:prepare(u) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then 
		print ("Newton solver failed at step "..step.."."); exit(); 
	end

	-- update new time
	time = solTimeSeries:time(0) + do_dt

	-- get oldest solution
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)

	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)

	RescaleDomainPointwise(dom)

	out:clear_selection()
	out:select(VelCmp, "velocity at begin" .. time)
	out:print("test", u,step,time)

	print(Integral(1.0, u, "left"));
	print(Integral(1.0, u, "bottom"));
end