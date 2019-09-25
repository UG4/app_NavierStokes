----------------------------------------------------------------------------------
--
--  Lua - Script to compute (stationary, 2d) the terminal sedimenting
--  velocity of a circular particle falling in stokes flow. It applies the
--  'MovingParticle' class as Navier Stokes discretisation implemented to
--  handle moving particle.
--  An analytical solution for the terminal velocity and the drag force
--  is provided by the literature for arbitrary radius of the particle,
--  see
--
--  Wagner, G. J., Moes, N., Liu, W. K., Belytschko, T.:
--  'The extended finite element method for rigid particles in Stokes flow'
--  Int. J. for Numerical Methods in Eng., (2001), pp. 293-313
--
--  Numerical results for the method implemented here are published in:
--
--  S. Höllbacher, G. Wittum:
--  'Rotational test spaces for a fully-implicit FVM and FEM for the DNS of
--  'fluid-particle interaction'
--  Journal of Computational Physics, vol. 393, (2019), pp. 186–213
--  open access link: https://authors.elsevier.com/sd/article/S0021999119303298
--  DOI: https://doi.org/10.1016/j.jcp.2019.05.004
--
--  and
--
--  S. Hoellbacher, G. Wittum:
--  'Gradient-consistent enrichment of finite element spaces for the DNS
--  'of fluid-particle interaction.'
--  submitted to Jounal of Computational Physics
--
--   Author: Susanne Hoellbacher
--
-----------------------------------------------------------------------------------


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Setup
--------------------------------------------------------------------------
--------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim 		 = util.GetParamNumber("-dim", 2, "world dimension")
numRefs 	 = util.GetParamNumber("-numRefs", 5, "number of grid refinements")
numPreRefs 	 = util.GetParamNumber("-numPreRefs", 5, "number of prerefinements (parallel)")

InitUG(dim, AlgebraType("CPU", dim+1));

if 	dim == 2 then
    gridName = util.GetParam("-grid", "grids/channel_4x12_irregular_tri_PrsBndCond.ugx")
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Parameter Setup
--------------------------------------------------------------------------
--------------------------------------------------------------------------

-- particle radius and center
radius = 0.1
center = {1.0, 0.0}

-- viscosity of the fluid
visc = 1.0

-- densities
densityFluid = 1.0
densityPrt = 2.0

-- gravitation constant (depending on the units used)
gravityConst = -9.81

-- boolian for choosing only stokes equation
stokes = true

-- for particle simulations, due to the irregular cutting of the particle interface
-- with the finite elements, the scheme is not symmetric. We need to set laplace = false
laplace = false

-- further settings to choose

volumeCompMode = true
-- default = true:  computes the volume analytically using the formluar for a circle or sphere;
-- false: computes the volume of the particle based on the cut elements, i.e. with polygonal boundary

StdFVAssembling = false
-- = true: uses non-projected, i.e. non-FlatTop ansatz- and test-spaces as in the original ficticious domain method

SkipNearInterfaceElem = 0
-- = 1: based on a given threshold elements are not considered as cut,
--      if their nodes are lying in a band of thickness equal the threshold around the interface
-- There are 2 possiblities to set a threshold:
-- movingParticle:set_threshold(lev, threshold):
--      --> sets the threshold for gridlevel = lev to the value 'threshold'
-- movingParticle:initialize_threshold(dom, numRefs, numRefs):
--      --> sets the to '0.25 * meanLength*meanLength', with meanLenght = MeanElementDiameter(domain, lev)


do_transfer = true
-- true: uses transfer operator adapted to the new coordinates of a cut element if
--       using a gmg as linear solver

Solver = 1
newtonMinDef =  util.GetParamNumber("-newtonMinDef", 1e-8)
newtonRed =  util.GetParamNumber("-newtonRed", 1e-9)

topLevel = numRefs
baseLevel = numRefs

--------------------------------------------------------------------------
--------------------------------------------------------------------------
--- Print Parameter
--------------------------------------------------------------------------
--------------------------------------------------------------------------

print(" Chosen Parameters:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("	   radius           = " .. radius)
print("    center           = " .. center[1].."  " .. center[2])
print("	   visc             = " .. visc)
print("	   densityFluid		= " .. densityFluid)
print("	   densityPrt		= " .. densityPrt)
if ( stokes ) then
print("    stokes           = true")
else
print("    stokes           = false")
end


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------
--------------------------------------------------------------------------

neededSubsets = {"Inner", "Inlet", "Outlet", "Upper", "Lower", "PressureBndCond"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Lagrange", 1) end
if dim >= 2 then approxSpace:add_fct("v", "Lagrange", 1) end

approxSpace:add_fct("p", "Lagrange", 1)

-- we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

OrderCuthillMcKee(approxSpace, true);

dbgWriter = GridFunctionDebugWriter(approxSpace)
--dbgWriter:set_vtk_output(true)

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- User Data Setup
--------------------------------------------------------------------------
--------------------------------------------------------------------------

function ZeroSolution(x, y, t)
    return 0.0
end


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Navier-Stokes Discretization Setup
--------------------------------------------------------------------------
--------------------------------------------------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
fctUsed = fctUsed .. ", p"

--------------------------------------------------------------------------
-- start loop over all radii
--------------------------------------------------------------------------
for i_rad = 0, 40, 1 do
radius = 0.1 + i_rad*0.01


NavierStokesDisc = NavierStokesFV1_cutElem(fctUsed, "Inner")

-- set upwind
fullUpwind = NavierStokesFullUpwind();
lpsUpwind = NavierStokesLinearProfileSkewedUpwind();

-- set stabilisation
stabFLOW = NavierStokesFLOWStabilization()
stabFIELDS = NavierStokesFIELDSStabilization()

stab = stabFIELDS
stab:set_upwind(fullUpwind)
stab:set_diffusion_length("RAW")

-- set up
NavierStokesDisc:set_stabilization(stab)
if not(stokes) then
    stab:set_upwind(lpsUpwind)
    NavierStokesDisc:set_upwind(lpsUpwind)
    NavierStokesDisc:set_pac_upwind(true)
end

NavierStokesDisc:set_peclet_blend(true)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_stokes(stokes)
NavierStokesDisc:set_laplace(laplace)

function Density(x, y, t)
    dens = densityFluid
    return dens
end

NavierStokesDisc:set_density("Density")
NavierStokesDisc:set_kinematic_viscosity(visc);


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Set up fluid-particle problem
--------------------------------------------------------------------------
--------------------------------------------------------------------------

-- set up domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)

--------------------------------------------------------------------------
-- Set up Particle functionality
--------------------------------------------------------------------------

-- add 1 particle with spherical shape to the discretisation
interfaceProvider = ParticleProviderSphere()
interfaceProvider:add(radius, MakeVec(center[1], center[2]), densityPrt)

cutElementHandler = CutElementHandler_FlatTop(dom:grid(), "u,v,p", interfaceProvider)

-- create the class handling the immersed particle
movingParticle = MovingParticle(domainDisc, NavierStokesDisc, cutElementHandler, densityFluid, visc)

movingParticle:set_gravity(true, gravityConst)
movingParticle:set_StdFV_assembling(StdFVAssembling)

if SkipNearInterfaceElem == 1 then
    -- sets the threshold to '0.25 * meanLength*meanLength', with meanLenght = MeanElementDiameter(domain, lev)
    movingParticle:initialize_threshold(dom, numRefs, numRefs)

    -- it is also possible, to set a user defined threshold = 'thres' by calling
    --   movingParticle:set_threshold(lev, thres)
end

interfaceProvider:print()

--------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------

--setup Walls
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("Upper,Lower, Outlet, Inlet")

--set pressure value in node
PressureBndCond = DirichletBoundary()
PressureBndCond:add("ZeroSolution", "p", "PressureBndCond")

-- set coupling condition on immersed boundary
InnerBndCond = movingParticle:get_BndCond()


-- add boundary conditions
domainDisc:add(InnerBndCond)
domainDisc:add(WallDisc)
domainDisc:add(PressureBndCond)


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Set up Solvers
--------------------------------------------------------------------------
--------------------------------------------------------------------------

-- create algebraic Preconditioner
ilu = ILU()
ilu:set_beta(-0.99)
--ilu:set_debug(dbgWriter)


baseSolver = BiCGStab()
baseSolver:set_preconditioner(ilu)
baseSolver:set_convergence_check(ConvCheck(10000000, 1e-11, 1e-11, false))
baseSolver:set_compute_fresh_defect_when_finished(true)

Base = LinearSolver()
Base:set_convergence_check(ConvCheck(2000000, 1e-13, 1e-13, true))
Base:set_preconditioner(ILU())

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(numPreRefs)
gmg:set_base_solver(SuperLU())
gmg:set_smoother(ilu)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
--gmg:set_debug(dbgWriter)

if do_transfer == true then
    transfer = ParticleTransfer(approxSpace, cutElementHandler)
    transfer:set_debug(dbgWriter)
    transfer:set_use_transposed(true)
    gmg:set_transfer(transfer)
end

-- create Linear Solver
convCheck = ConvCheck(1000, 1e-10, 1e-10, true)
convCheck:set_verbose(true)

linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)
linSolver:set_compute_fresh_defect_when_finished(true)
--linSolver:set_debug(dbgWriter)

-- create Exact solver
exactSolver = SuperLU()

-- choose a solver
if Solver == 1 then solver = exactSolver
else                solver = linSolver
end

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(1000000)
newtonConvCheck:set_minimum_defect(newtonMinDef)
newtonConvCheck:set_reduction(newtonRed)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(10)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
--newtonSolver:set_debug(dbgWriter)


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------
--------------------------------------------------------------------------

op = AssembledOperator(domainDisc)
op:init()

u = GridFunction(approxSpace)
u:set(0)

function Pressure_StartValue2d(x, y, t) return 0.0 end
function VelX_StartValue2d(x, y, t)     return 0.0 end
function VelY_StartValue2d(x, y, t)     return 0.0 end

time = 0.0
Interpolate("Pressure_StartValue"..dim.."d", u, "p", time);
Interpolate("VelX_StartValue"..dim.."d", u, "u", time);
Interpolate("VelY_StartValue"..dim.."d", u, "v", time);

-- call adjust_solution() to set the soltuion in the particle domain
domainDisc:adjust_solution(u)


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Init, prepare, apply newton Solver
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

-- we have to prepare the particle class. This sets the marker for the
-- cut elements and initializes the values
-- initialize 'movingIParticle' class with topLevel as argument TWICE:
--  ==> no loop for multigrid
movingParticle:init(u, approxSpace, baseLevel, topLevel)


-- init and prepare newton
newtonSolver:init(op)
if newtonSolver:prepare(u) == false then
    print ("Newton solver prepare failed."); exit();
end

-- apply solver
if newtonSolver:apply(u) == false then
    print ("Newton solver apply failed."); exit();
end



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Output stuff
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

domainDisc:adjust_solution(u)


-- write particle velocity to file
filename_output_vel = "./SedimentingStationary_all_radii_numRefs"..numRefs..".txt"
movingParticle:print_velocity(u, topLevel, time, filename_output_vel)


print("done computation for radius = "..radius.."\n")

--------------------------------------------------------------------------------

-- end of radius-loop 'i_rad'
end

print("done.")


