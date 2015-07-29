--------------------------------------------------------------------------------
--
--  Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--  Author: Josef Dubsky, Andreas Vogel
--
--------------------------------------------------------------------------------

-- Execute, e.g., via `mpirun -n 4 bin/ugshell -ex apps/navier_stokes/cylinder/main.lua`

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("../navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

dim             = util.GetParamNumber("--dim", 2, "world dimension")
numRefs	        = util.GetParamNumber("--numRefs", 2, "number of grid refinements")
numPreRefs      = util.GetParamNumber("--numPreRefs", 0, "number of prerefinements (parallel)")
bConvRates      = util.HasParamOption("--convRate", "compute convergence rates")
bBenchmarkRates = util.HasParamOption("--benchRate", "compute benchmark rates")

order           = util.GetParamNumber("--order", 1, "order pressure and velocity space")
vorder          = util.GetParamNumber("--vorder", order, "order velocity space")
porder          = util.GetParamNumber("--porder", order-1, "order pressure space")

bStokes         = util.HasParamOption("--stokes", "If defined, only Stokes Eq. computed")
bNoLaplace      = util.HasParamOption("--nolaplace", "If defined, only laplace term used")
bExactJac       = util.HasParamOption("--exactjac", "If defined, exact jacobian used")
bPecletBlend    = util.HasParamOption("--pecletblend", "If defined, Peclet Blend used")
upwind          = util.GetParam("--upwind", "lps", "Upwind type")
stab            = util.GetParam("--stab", "flow", "Stabilization type")
diffLength      = util.GetParam("--difflength", "cor", "Diffusion length type")

discType, vorder, porder = util.ns.parseParams()

Viscosity = 1e-3
Um = 0.3
if dim == 3 then Um = 0.45 end
H = 0.41
L = 0.1
Umean2 = math.pow(2/3*Um, 2)

ref = {}
ref.CD = 5.57953523384
ref.CL = 0.010618948146
ref.DeltaP = 0.11752016697

if dim == 2 then 
	gridName = util.GetParam("--grid", "../grids/cylinder.ugx")
	--gridName = util.GetParam("--grid", "grids/box.ugx")
	--gridName = util.GetParam("--grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("--grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("--grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("--grid", "grids/cylinder_rotate_box_tri_fine.ugx")
elseif dim == 3 then
	gridName = util.GetParam("--grid", "grids/cylinder3d.ugx")
    --gridName = util.GetParam("--grid", "grids/cylinder3d_fine.ugx")
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

ug_load_script("domain.lua")

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

ug_load_script("discretization.lua")

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

ug_load_script("solution.lua")

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

ug_load_script("run.lua")
