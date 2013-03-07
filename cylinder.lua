--------------------------------------------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--   Author: Josef Dubsky, Andreas Vogel
--
--------------------------------------------------------------------------------

-- Right at the beginning we load a lot of util functions, that help us basically
-- to programm a domain independent lua-script and provide some basic helper
-- functions. Most the util functions are in the namespace util, i.e. they
-- are used by 'util.SomeFunction'
ug_load_script("ug_util.lua")

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here. The dimesion is
-- read in from the bash command line passing e.g. the option "-dim 2". If no
-- dim-option is specified the default value in the second argument of
-- GetParamNumber is used.
dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- Next, we have to choose some algebra. ug4 provides several algebra, e.g.
-- there are block structured matrices or simple double-valued matrices. We
-- decide to use the double-valued CSR Matrix. This is the default case for the
-- Algebra chooser and so we leave the intiallizer of the AlgebraChooser empty.
InitUG(dim, AlgebraType("CPU", 1));

-- Next, we decide which grid to use. This can again be passed as a command line
-- option or a default value is used.
if 	dim == 2 then gridName = util.GetParam("-grid", "grids/cylinder.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numTotalRefs = util.GetParamNumber("-numRefs", 0)

-- Calculate the number of post-refs and make sure that the result makes sense.
numPostRefs = numTotalRefs - numPreRefs
if numPostRefs < 0 then
	print("WARNING:\tnumPreRefs exceeds the number of total refinements.")
	print("\t\t\tNo refinement will be preformed after distribution.")
	numPostRefs = 0
end

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numTotalRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Create the domain and load a grid
dom = Domain()
LoadDomain(dom, gridName)

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)

-- We will create a refiner now. This is a tool associated with a domain.
-- UG defines factory methods for refiners, which automatically choose
-- the right refiner for the given context, i.e. different refiners are
-- created depending on whether we are in a parallel or in a serial environment.
-- Note that another factory method is HangingNodeDomainRefiner, which is
-- subject to a later tutorial.
refiner = GlobalDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
	refiner:refine()
end

-- Distribute the refined domain to all involved processes
if util.DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end

-- perform post-refinement
for i = 1, numPostRefs do
	refiner:refine()
end

--------------------------------------------------------------------------------
-- Approximation Space
--------------------------------------------------------------------------------

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

-- Lets define a list of all subsets that we need
neededSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall", "CylinderWall"}

-- Now we loop all subsets an search for it in the SubsetHandler of the domain
if util.CheckSubsets(dom, neededSubsets) == false then print("Wrong subsets detected.") end

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"}; FctCmp = {"u", "p"};
elseif  dim == 2 then VelCmp = {"u", "v"}; FctCmp = {"u", "v", "p"};
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"};
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- we add the velocity and pressure as Lagrange Ansatz function of first order
approxSpace:add_fct(FctCmp, "Lagrange", 1) 

-- finally we print some statistic on the distributed dofs
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

-- We set up the discretization. The basic idea is to first create the single
-- components of the discretization (e.g. element discretization, boundary
-- conditions, ...) and then glue everything together in one DomainDiscretization
-- object that does nothing else than grouping the single components. When it
-- comes to assembling this object runs all scheduled discretization component
-- one after the other.

-- Now we create the Navier-Stokes disc. Here we use a util function, that
-- creates the disc for the right world dimension according to the loaded domain
-- (this is handled internally, since the approxSpace knows the domain and the
-- domain knows the dimension). We tell the elem Disc to assemble on all elements
-- that are in the subset "Inner". If we would have several domains, where we
-- would like to do the same, this could be done by passing a list of subsets
-- separated by ',', (e.g. "Inner1, Inner2, Inner3").
NavierStokesDisc = NavierStokes(FctCmp, {"Inner"})

-- Now, we have to setup the stabilization, that is used for the Continuity Equation.
-- The stabilization is passed to the Navier-Stokes elem disc as an object.
-- Therefore, we created it now and will pass it to the disc. But first, we have
-- to create the Upwind scheme, that is used inside the stabilization. There are
-- several possibilities:

noUpwind = NavierStokesNoUpwind();
fullUpwind = NavierStokesFullUpwind();
skewedUpwind = NavierStokesSkewedUpwind();
LPSUpwind = NavierStokesLinearProfileSkewedUpwind();
POSUpwind = NavierStokesPositiveUpwind();

-- Now, we create the stabilization ...
fieldsStab = NavierStokesFIELDSStabilization()

-- ... and set the upwind
fieldsStab:set_upwind(fullUpwind)

-- We also can choose, how the diffusion length of the stabilization is computed.
-- Under the option we pick on:
fieldsStab:set_diffusion_length("RAW")
--fieldsStab:set_diffusion_length("FIVEPOINT")
--fieldsStab:set_diffusion_length("COR")

-- Next we set the options for the Navier-Stokes elem disc ...
NavierStokesDisc:set_stabilization(fieldsStab)
NavierStokesDisc:set_conv_upwind(fullUpwind)
NavierStokesDisc:set_peclet_blend(false)
NavierStokesDisc:set_exact_jacobian(false)

-- ... and finally we choose a value for the kinematic viscosity.
NavierStokesDisc:set_kinematic_viscosity(1.0e-1);


-- Next, lets create the boundary conditions. We will use lua-callback functions
-- to implement the values of the dirichlet-bnd conditions. So lets define
-- some functions in lua: Parameters are the coordinates of the point at which
-- the value should be evaluated and the time of the current time-step.
-- The dirichlet boundary callback has to return two values: a boolean
-- defining whether the given point really is a dirichlet boundary point
-- and the boundary value itself.
-- Note that we use math-functions from luas standard math-library.
-- (here the . operator is used, since math is not an object but a library)
function inletVelX2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return true, 4 * Um * y * (H-y) / (H*H)
end

function inletVelY2d(x, y, t)
	local H = 0.41
	local Um = 0.2
--	return true, Um * math.sin(y/H*2*3.1415)
	return true, 0.0
end

-- We need neumann boundary for the pressure at Inlet
function inletPressure2d(x, y, t)
	local b, vel = inletVelX2d(x,y,t)
	return true, -1.0 * vel
end

function inletVel2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return 4 * Um * y * (H-y) / (H*H), 0.0
end

-- setup outlet
-- use "no normal stress" outflow condition
OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
OutletDisc:add("Outlet")

-- alternatively use "zero pressure" outflow condition
-- OutletDisc = DirichletBoundary()
-- OutletDisc:add(0.0, "p", "Outlet")

-- setup inlet
InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Inlet")

-- setup wall
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("UpperWall,LowerWall,CylinderWall")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)
domainDisc:add(WallDisc)
domainDisc:add(OutletDisc)

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

-- Now lets solve the problem. Create a vector of unknowns and a vector
-- which contains the right hand side. We will use the approximation space
-- to create the vectors. Make sure to create the vector for the same
-- dofs as set to linOp through linOp:set_dof_distribution.
-- Note that the vectors automatically have the right size.
u = GridFunction(approxSpace)

-- Init the vector representing the unknowns with 0 to obtain
-- reproducable results.
u:set(0)
time = 0.0

-- we need a linear solver that solves the linearized problem inside of the
-- newton solver iteration. We create an exact LU solver here and an HLibSolver.
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(LU())
gmg:set_smoother(ILU())
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)

-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(1000, 1e-7, 1e-3, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(1000, 1e-7, 1e-3, true))

-- choose a solver
solver = gmgSolver

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver(domainDisc)
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(ConvCheck(20, 1e-6, 1e-10, true))
newtonSolver:set_line_search(StandardLineSearch(20, 1.0, 0.5, true))
newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- Now we can apply the newton solver. A newton itertation is performed to find
-- the solution.
if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

-- Finally we're nearly done. The only thing left to do is to write the
-- solution to a file which can then be examined using e.g. Paraview.
-- (Open "Solution.vtu" in paraview to view the complete domain
vtkWriter = VTKOutput()
vtkWriter:select_all(false)
vtkWriter:select_nodal(VelCmp, "velocity")
vtkWriter:select_nodal("p", "p")
vtkWriter:print("Solution", u)