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

local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/?.lua"

ug_load_script("ug_util.lua")
-- ug_load_script("util/domain_disc_util.lua")
-- ug_load_script("util/conv_rates_static.lua")
ug_load_script("plugins/Limex/limex_util.lua")
ug_load_script("navier_stokes_util.lua")

RequiredPlugins({"Limex", "NavierStokes"})




-- Command line arguments
local ARGS = {
  dim     = util.GetParamNumber("-dim", 2, "world dimension"),
  doSteadyState = util.HasParamOption("--steady-state", "Compute steady state solution"),
  bStokes   = util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed"),
  bExactJac   = util.HasParamOption("-exactjac", "If defined, exact jacobian used"),
  bLaplace  = not util.HasParamOption("-nolaplace", "If defined, only laplace term used"),
  bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used"),
  
  upwind      = util.GetParam("-upwind", "lps", "Upwind type"),
  
  stab        = util.GetParam("-stab", "flow", "Stabilization type"),
  diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type"),
  
  stabGrad             = util.GetParamNumber("--stabGrad", 0.1, "Stabilization parameter."),
  stabStreamline      = util.GetParamNumber("--stabStreamline", 0.0, "Stabilization parameter."),
  stabDiv             = util.GetParamNumber("--stabDiv", 0.0, "Stabilization parameter."),
 

  limexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages: 0:steady state, 1:std. integration, <= 2:LIMEX"),
  limexDebugLevel = util.GetParamNumber("--limex-debug-level", 2, "debug level"),
  limexTOL = util.GetParamNumber("--limex-tol", 1e-2, "debug level"),
  
  solverID =  util.GetParam("--solver-id", "gmg", "superlu, gmg"),
  
  numRefs   = util.GetParamNumber("--numRefs", 1, "number of grid refinements"),
  numPreRefs  = util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
}

ARGS.discType, ARGS.vorder, ARGS.porder = util.ns.parseParams()





-- Shortcut.
local dim =ARGS.dim

if dim == 2 then 
	ARGS.gridName = util.GetParam("-grid", "grids/cylinderg.ugx")
	--gridName = util.GetParam("-grid", "grids/box.ugx")
	--gridName = util.GetParam("-grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_rotate_box_tri_fine.ugx")
elseif dim == 3 then
	ARGS.gridName = util.GetParam("-grid", "grids/cylinder3d.ugx")
--	gridName = util.GetParam("-grid", "grids/cylinder3d_fine.ugx")
else print("Selected Dimension not supported. Exiting."); exit(); end


-- Lets write some info about the choosen parameter
print(" Selected Parameter:")
print("    dim              = " .. dim)

print("    numTotalRefs     = " .. ARGS.numRefs)
print("    numPreRefs       = " .. ARGS.numPreRefs)
print("    grid             = " .. ARGS.gridName)

print("    porder           = " .. ARGS.porder)
print("    vorder           = " .. ARGS.vorder)
print("    type             = " .. ARGS.discType)
print("    only stokes      = " .. tostring(ARGS.bStokes))
print("    only laplace     = " .. tostring(ARGS.bLaplace))
print("    exact jacobian   = " .. tostring(ARGS.bExactJac))
print("    peclet blend     = " .. tostring(ARGS.bPecletBlend))
print("    upwind           = " .. ARGS.upwind)
print("    stab             = " .. ARGS.stab)
print("    stabGrad         = " .. ARGS.stabGrad)
print("    stabDiv          = " .. ARGS.stabDiv)
print("    stabStreamline   = " .. ARGS.stabStreamline)
print("    diffLength       = " .. ARGS.diffLength)

print("    LIMEX.TOL        = " .. ARGS.limexTOL)
print("    LIMEX.NumStages  = " .. ARGS.limexNStages)

--------------------------------------------------------------------------------
-- Debug output.
--------------------------------------------------------------------------------
local logAssistant = GetLogAssistant()
logAssistant:set_debug_level("LIB_LIMEX", ARGS.limexDebugLevel)
--logAssistant:set_debug_levels(10)
--logAssistant:set_debug_level("LIB_DISC_MULTIGRID", 10)

--------------------------------------------------------------------------------
-- InitUG
--------------------------------------------------------------------------------
InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Problem setup.
--------------------------------------------------------------------------------
local myProblem=require("cylinder-limex-config")
myProblem:Init(ARGS) 

local dom = myProblem:CreateDomain(ARGS.gridName, ARGS.numRefs, ARGS.numPreRefs)
local approxSpace = myProblem:CreateApproxSpace(dom)
local domainDisc = myProblem:CreateDomainDisc(approxSpace)
-- local solver = myProblem:CreateSolver(approxSpace)

local myVelCmp = myProblem:GetVelocityCmps()
local myDbgWriter = GridFunctionDebugWriter(approxSpace)
-- myDbgWriter.set_conn_viewer_output(true)

-- local transfer = 
-- local transfer = StdTransfer()
-- transfer:enable_p1_lagrange_optimization(true)
-- transfer:set_debug(myDbgWriter)

--------------------------------------------------------------------------------
-- Solver setup.
--------------------------------------------------------------------------------
local sscSmootherDesc = { type = "ssc", vertex = {{"p"}, {"u","v"}} }
local mySmootherDesc = sscSmootherDesc

local mySolverDesc = {

    type = "newton",
    lineSearch = StandardLineSearch(10, 1.0, 0.5, true, true),
    convCheck = ConvCheck(50, 1e-13, 1e-8, true, false),-- "standard",  -- verbose, suppress

    linSolver = {
      type = "bicgstab", --"linear",
      convCheck = ConvCheck(100, 1e-16, 1e-8, true),
    }
}

-- Special part (based on command line args).
if (ARGS.solverID == "superlu") then
  mySolverDesc.linSolver.type = "superlu"
elseif (ARGS.solverID == "gmg") then

  mySolverDesc.linSolver.type = "bicgstab"
  mySolverDesc.linSolver.precond =  {
    type    = "gmg",

    approxSpace = approxSpace,
    baseLevel = 0,
    baseSolver  = "superlu",
    cycle = "F", 
    rap = true,
    -- discretization = domainDisc,

    preSmooth = 1,
    postSmooth = 1,
    smoother  =  mySmootherDesc,

    transfer = "std",

    -- transfer = transfer, 
    debug = true,
    
    -- [[ 
    debugSolver = {

      type = "linear",
      convCheck = ConvCheck(20, 5e-12, 1e-99, true),
      precond =  mySmootherDesc,

      approxSpace = approxSpace,
    }, -- debugSolver
  --]]
  }  


else  

  mySolverDesc.linSolver.precond = mySmootherDesc
  --{ type = "ssc", damping = 1.0, vertex = {{"p"}, myVelCmp} } --]]

end


local solver = util.solver.CreateSolver(mySolverDesc, solverutil)
--solver:set_debug(GridFunctionDebugWriter(approxSpace))
print(solver:config_string())


--------------------------------------------------------------------------------
-- LIMEX configuration
--------------------------------------------------------------------------------
local myLimexDesc = {}
    myLimexDesc["STAGES"] = ARGS.limexNStages -- tolerance
    myLimexDesc["TOLERANCE"] = ARGS.limexTOL-- tolerance
    myLimexDesc["REDUCTION"] = 0.5  -- reduction
    myLimexDesc["INCREASE"]  = 2.0 -- increase of time step
    myLimexDesc["SAFETY"] = 0.8    -- safety factor
    -- local errorEst = myLimexDesc["ESTIMATOR"]

    myLimexDesc["SPACES"] = {
      H1SemiComponentSpace("u", 4), 
      H1SemiComponentSpace("v", 4), 
    --  L2QuotientSpace("p", 2), -- need to factor out constants!
    }

    -- myLimexDesc["DEBUG"] = myDbgWriter -- enable output of matrices

--------------------------------------------------------------------------------
-- VTK writer.
--------------------------------------------------------------------------------
local vtkWriter = VTKOutput()
vtkWriter:select("p", "pressure")
vtkWriter:select(myVelCmp, "velocity")

--------------------------------------------------------------------------------
-- Solve problem.
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
myProblem:SetInitialValues(u)

if (ARGS.doSteadyState) then
  -- Steady state solution.
  myProblem:ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)
  vtkWriter:print("CylinderSteadyState", u)
else
  -- Transient solution.
  local cTransient = 
    {
      -- Start and stop time --
      tStart = 0.0,
      tStop = 8.0,

      -- Fractional-step-theta (same amount of work as CN, but stable) --
      --  scheme = "fracstep", maxStepSize = 0.1, minStepSize = 8.0/4096.0, redStepSize = 0.5,

      -- Limex time stepping
      scheme = "limex", maxStepSize = 0.5, minStepSize = (8.0/2048.0)*1e-8, redStepSize = 0.5, 
     
      -- Implicit Euler (requires a small time step) --
      --scheme = "impleuler", maxStepSize = 0.001, minStepSize = 0.0001, redStepSize = 0.5

      -- SDIRK (3rd-order) a.k.a. Alexander3 (we call it `sdirk3`) --
      --scheme = "sdirk", orderOrTheta = 3, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 

      -- SDIRK (4th-order) a.k.a. Hairer, Wanner, L-stable DIRK (we call it `sdirk4`), TODO: Check --
      --scheme = "sdirk", orderOrTheta = 4, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 
    }



  local ComputeEffectiveQuantities = function (u, step, time) myProblem:EvalIntegralQuantities2D(u, step, time) end

  if (cTransient.scheme == "limex") then 


    util.SolveNonlinearProblemLimex(
      u, domainDisc, solver,
      vtkWriter, "limex_solution",
      cTransient.tStart, 
      cTransient.tStop, 
      cTransient.maxStepSize/10.0,  -- start step size
      cTransient.minStepSize, 
      cTransient.maxStepSize, 
      myLimexDesc, 
      ComputeEffectiveQuantities)
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

