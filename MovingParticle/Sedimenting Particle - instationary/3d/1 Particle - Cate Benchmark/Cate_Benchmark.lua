
-----------------------------------------------------------------------------------
--
--  Lua - Script to compute a sedimenting 3d-particle with spherical shape. It
--  applies the immersed interface method provided by the 'MovingInterface' class.
--  Experimental and numerical reference data are provided by:
--
--  A. ten Cate, C. H. Nieuwstad, J. J. Derksen, and H. E. A. Akker.
--  'Particle imaging velocimetry experiments and lattice-Boltzmann simulations'
--  'on a single sphere settling under gravity.'
--  Physics of Fluids, vol. 14, 2002, pp. 4012-4025
--
--
--  Numerical results for the method implemented here are published in:
--
--  S. Hoellbacher, G. Wittum:
--  'Gradient-consistent enrichment of finite element spaces for the DNS
--  'of fluid-particle interaction.'
--  submitted to Jounal of Computational Physics
--
--
--   Author: Susanne Hoellbacher
--
-----------------------------------------------------------------------------------


-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-- Setup
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

dim 		 = util.GetParamNumber("-dim", 3, "world dimension")
numRefs 	 = util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	 = util.GetParamNumber("-numPreRefs", 2, "number of prerefinements (parallel)")
dt           = util.GetParamNumber("-dt", 0.004)
numTimeSteps = util.GetParamNumber("-numTimeSteps", 1250)   -- T = 5.0

InitUG(dim, AlgebraType("CPU", dim+1));

if 	dim == 3 then
    gridName = util.GetParam("-grid", "Sedimenting_10x10x16_regular.ugx")
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end


-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-- Parameter Setup
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

-- particle radius and center
radius = 0.075
center = {0.0, 0.0, 1.2}

-- densities
densityFluid = 0.970
densityPrt = 1.12

-- viscosity of the fluid
visc_mu = 0.0373
visc = visc_mu/densityFluid

-- gravitation constant (depending on the units used)
gravityConst = -98.1

-- boolian for choosing only stokes equation
stokes = false

-- for particle simulations, due to non-symmetric sub-control-volume-faces along
-- the interface, the scheme is not symmetric, so we need to set laplace = false
laplace = false

-- further settings to choose

StdFVAssembling = false
-- = true: uses non-projected, i.e. non-Flat-Top ansatz- and test-spaces

print_cutElemData = true
-- = true: prints computed coordinates of the corners of cut elements and their local indices to a file

SkipNearInterfaceElem = 0
-- = 1: based on a given threshold elements are not considered as cut,
--      if their nodes are lying in a band of thickness equal the threshold around the interface
-- There are 2 possiblities to set a threshold:
-- movingParticle:set_threshold(lev, threshold):
--      --> sets the threshold for gridlevel = lev to the value 'threshold'
-- movingParticle:initialize_threshold(dom, numRefs, numRefs):
--      --> sets the to '0.25 * meanLength*meanLength', with meanLenght = MeanElementDiameter(domain, lev)

do_transfer = false
-- true: uses transfer operator adapted to the new coordinates of a cut element if
--       using a gmg as linear solver

Solver = 1
newtonRed =  util.GetParamNumber("-newtonRed", 1e-10)

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
print("    dt               = " .. dt)
print("    numTimeSteps     = " .. numTimeSteps)
print("    grid             = " .. gridName)
print("	   radius           = " .. radius)
print("    center           = " .. center[1].."  " .. center[2].."  " .. center[3])
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

dom = Domain()
LoadDomain(dom, gridName)

approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Lagrange", 1) end
if dim >= 2 then approxSpace:add_fct("v", "Lagrange", 1) end
if dim >= 3 then approxSpace:add_fct("w", "Lagrange", 1) end

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
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokesFV1_cutElem(fctUsed, "Inner")

-- set upwind
fullUpwind = NavierStokesFullUpwind();
lpsUpwind = NavierStokesLinearProfileSkewedUpwind();

-- set stabilisation
stabFLOW = NavierStokesFLOWStabilization()
stabFIELDS = NavierStokesFIELDSStabilization()

stabFLOW:set_upwind(lpsUpwind)
stabFLOW:set_diffusion_length("RAW")

-- set up
NavierStokesDisc:set_upwind(lpsUpwind)
NavierStokesDisc:set_stabilization(stabFLOW)
NavierStokesDisc:set_pac_upwind(true)
NavierStokesDisc:set_peclet_blend(false)
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
interfaceProvider:add(radius, MakeVec(center[1], center[2], center[3]), densityPrt)

cutElementHandler = CutElementHandler_FlatTop(dom:grid(), "u,v,w,p", interfaceProvider)

-- create the class handling the immersed particle
movingParticle = MovingParticle(domainDisc, NavierStokesDisc, cutElementHandler, densityFluid, visc)

movingParticle:set_gravity(true, gravityConst)
movingParticle:set_time_step(dt)
movingParticle:set_StdFV_assembling(StdFVAssembling)

if SkipNearInterfaceElem == 1 then
    -- sets the threshold to '0.25 * meanLength*meanLength', with meanLenght = MeanElementDiameter(domain, lev)
    movingParticle:initialize_threshold(dom, numRefs, numRefs)

    -- it is also possible, to set a user defined threshold 'thres' by calling
    --   movingParticle:set_threshold(lev, thres)
end

-- prints computed coordinates of all corners of cut elements and their local indices to a file
-- see './CutElementData.txt'
movingParticle:set_print_cutElemData(print_cutElemData)

-- print particle data to terminal
interfaceProvider:print()


local ref = GlobalDomainRefiner(dom)

if NumProcs() > 1 then
	pu = ParticleUnificator(dom)
	cdgm = ClusteredDualGraphManager()
	cdgm:add_unificator(SiblingUnificator())
	cdgm:add_unificator(pu)

	balancer.partitioner = "parmetis"
	balancer.staticProcHierarchy = true
	balancer.firstDistLvl = 1
	balancer.firstDistProcs = 16
	balancer.redistProcs = 4
	balancer.redistSteps = 1
	balancer.parallelElementThreshold = 4
	balancer.qualityThreshold = 1.1
	balancer.ParseParameters()
	bal = balancer.CreateLoadBalancer(dom, balancerDesc)
	balancer.defaultPartitioner:set_dual_graph_manager(cdgm)

	pu:update_particles(interfaceProvider)
	  
	for i = 1, numPreRefs do
		ref:refine()
		bal:rebalance()
	end


	for i = numPreRefs+1, numRefs do
		ref:refine()
		bal.rebalance("util.refinement: adaption-" ..i)
	end

	if verbose and NumProcs() > 1 then
		bal.print_quality_records()
	end

else
	for i = 1, numPreRefs do
		ref:refine()
	end
	for i = numPreRefs+1, numRefs do
		ref:refine()
	end
end

movingParticle:set_element_diameter(MaxElementDiameter(dom, topLevel))

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------

--setup Walls
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("Front,Back,Top,Bottom,Left,Right,PressureBndCond")

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
-- Solution of the Problem
--------------------------------------------------------------------------
--------------------------------------------------------------------------

time = 0.0

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)                           -- 1.0 is implicit euler
op = AssembledOperator(timeDisc)

op:init()

u = GridFunction(approxSpace)
u:set(0)


-- call adjust_solution() to set the soltuion in the particle domain
domainDisc:adjust_solution(u)

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Set up Solvers
--------------------------------------------------------------------------
--------------------------------------------------------------------------

-- create algebraic Preconditioner
ilu = ILU()
--ilu:set_debug(dbgWriter)

baseSolver = BiCGStab()
baseSolver:set_preconditioner(ilu)
baseSolver:set_convergence_check(ConvCheck(10000000, 1e-11, 1e-11, false))
baseSolver:set_compute_fresh_defect_when_finished(true)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(numPreRefs)
gmg:set_base_solver(baseSolver)
gmg:set_smoother(ilu)
--gmg:set_smoother(Jacobi(0.8))  -- Jacobi schlechter fuer FT!
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
--exactSolver:set_minimum_for_sparse(10000000)

-- choose a solver
if NumProcs() > 1 then
	solver = LinearSolver()
	solver:set_preconditioner(baseSolver)
	solver:set_convergence_check(convCheck)
	solver:set_compute_fresh_defect_when_finished(true)
	--solver:set_debug(dbgWriter)
else
	if Solver == 1 then solver = exactSolver
	else                solver = linSolver
	end
end

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(2000)
newtonConvCheck:set_minimum_defect(1e-8)
newtonConvCheck:set_reduction(newtonRed)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(10)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.9)
newtonLineSearch:set_accept_best(true)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
--newtonSolver:set_debug(dbgWriter)


--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Init, prepare newton Solver
--------------------------------------------------------------------------
--------------------------------------------------------------------------


newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
    print ("Newton solver prepare failed."); exit();
end

-- we have to prepare the particle class. This sets the marker for the
-- cut elements and initializes the values
-- initialize 'movingIParticle' class with topLevel as argument TWICE:
--  ==> no loop for multigrid
movingParticle:init(u, approxSpace, topLevel, topLevel)

prtIndex = 0
numCutElem = movingParticle:get_numCutElements(topLevel, prtIndex)
print("\n #cut elements: " .. numCutElem)
print("\n For detailed info on the cut element geometry see file ./'CutElementData.txt'\n")


SaveVectorForConnectionViewer(u, "StartSolution.vec")

--------------------------------------------------------------------------
--------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------
--------------------------------------------------------------------------

-- start
time = 0.0
step = 0


vtkWriter = VTKOutput()
vtkWriter:select_all(false)
vtkWriter:select_nodal("u,v,w", "velocity")
vtkWriter:select_nodal("p", "pressure")
vtkWriter:print("Sedimenting_Sphere_GL"..numRefs, u, step, time, false)

-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


--------------------------------------------------------------------------
-- start time iteration
--------------------------------------------------------------------------

for step = 1, numTimeSteps do

    print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

    do_dt = dt

    timeDisc:prepare_step(solTimeSeries, do_dt)

    -- prepare newton solver
    if newtonSolver:prepare(u) == false then
        print ("Newton solver failed at step "..step.."."); exit();
    end


    -- apply newton solver
    if newtonSolver:apply(u) == false then
        print ("Newton solver failed at step "..step.."."); exit();
    end

    filename_output_vel = "./Sedimenting_Sphere_GL"..numRefs..".txt"
    movingParticle:print_velocity(u, topLevel, time, filename_output_vel)


    -- update new time
    time = solTimeSeries:time(0) + do_dt
    --vtkWriter:print("Sedimenting_Sphere_GL"..numRefs, u, step, time, false)

    interfaceProvider:print()

--------------------------------------------------------------------------
-- update particle data
--------------------------------------------------------------------------

-- updates data on cut elements; updates the coordinates of the paricle;
-- stores particle solution of last time step for access during assembling
	if NumProcs() > 1 then
		movingParticle:pre_balancing_update(u, approxSpace, baseLevel, topLevel, time, do_dt)
		pu:update_particles(interfaceProvider)     
		bal:rebalance()
		movingParticle:post_balancing_update(u, approxSpace, baseLevel, topLevel, time, do_dt)
	else
	    movingParticle:update(u, approxSpace, baseLevel, topLevel, time, do_dt)
	end


-- prints particle info to terminal
    interfaceProvider:print()


--------------------------------------------------------------------------
-- get oldest solution
--------------------------------------------------------------------------

    oldestSol = solTimeSeries:oldest()

    VecScaleAssign(oldestSol, 1.0, u)
    solTimeSeries:push_discard_oldest(oldestSol, time)

    print("++++++ TIMESTEP " .. step .. "  END ++++++");


end -- end time-loop

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- end time iteration
--------------------------------------------------------------------------
--------------------------------------------------------------------------


tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.");

domainDisc:adjust_solution(u)

SaveVectorForConnectionViewer(u, "FinalSolution.vec")

vtkWriter = VTKOutput()
vtkWriter:select_all(false)
vtkWriter:select_nodal("u,v,w", "velocity")
vtkWriter:select_nodal("p", "pressure")
vtkWriter:print("Sedimenting_Disc_GL"..numRefs, u)


print("done.")
