------------------------------------------------------------------
--
--   Lua - Script to test the Navier-Stokes implementation
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and is intended to test the implementation.
--
--   Author: Raphael Prohl, Andreas Vogel
--
-------------------------------------------------------------------

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
InitUG(dim, AlgebraType("CPU", dim+1));

-- Next, we decide which grid to use. This can again be passed as a command line
-- option or a default value is used.
if 	dim == 2 then gridName = util.GetParam("-grid", "grids/channel20Q18x10.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs", 0)

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Create the domain and load a grid
neededSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

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
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

OrderLex(approxSpace, "lr");
--OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

-- We set up the discretization. The basic idea is to first create the single
-- components of the discretization (e.g. element discretization, boundary
-- conditions, ...) and then glue everything together in one DomainDiscretization
-- object that does nothing else than grouping the single components. When it
-- comes to assembling this object runs all scheduled discretization component
-- one after the other.

-- Lets begin with the element discretization. This is the core part for the
-- Navier-Stokes equation and provides the local stiffness- and mass matrices,
-- that are then added to the global matrices. 

-- Now we create the Navier-Stokes disc. Here we use a util function, that
-- creates the disc for the right world dimension according to the loaded domain
-- (this is handled internally, since the approxSpace knows the domain and the
-- domain knows the dimension). We tell the elem Disc to assemble on all elements
-- that are in the subset "Inner". If we would have several domains, where we
-- would like to do the same, this could be done by passing a list of subsets
-- separated by ',', (e.g. "Inner1, Inner2, Inner3").
elemDisc = NavierStokes(FctCmp, {"Inner"})

-- Now, we have to setup the stabilization, that is used for the Continuity Equation.
-- The stabilization is passed to the Navier-Stokes elem disc as an object.
-- Therefore, we created it now and will pass it to the disc. But first, we have
-- to create the Upwind scheme, that is used inside the stabilization. There are
-- several possibilities:

--upwind = NavierStokesNoUpwind();
--upwind = NavierStokesFullUpwind();
--upwind = NavierStokesSkewedUpwind();
upwind = NavierStokesLinearProfileSkewedUpwind();
--upwind = NavierStokesRegularUpwind();
--upwind = NavierStokesPositiveUpwind();

-- Now, we create the stabilization ...
stab = NavierStokesFIELDSStabilization()
--stab = NavierStokesFLOWStabilization()

-- ... and set the upwind
stab:set_upwind(upwind)

-- We also can choose, how the diffusion length of the stabilization is computed.
-- Under the option we pick on:
--stab:set_diffusion_length("RAW")
stab:set_diffusion_length("FIVEPOINT")
--stab:set_diffusion_length("COR")

-- Next we set the options for the Navier-Stokes elem disc ...
elemDisc:set_stabilization(stab)
elemDisc:set_conv_upwind(upwind)
elemDisc:set_peclet_blend(false)
elemDisc:set_exact_jacobian(false)
elemDisc:set_stokes(true)
elemDisc:set_laplace(true)

-- ... and finally we choose a value for the kinematic viscosity.
elemDisc:set_kinematic_viscosity(0.01);

-- Next, lets create the boundary conditions. Lets define
-- some functions in lua: Parameters are the coordinates of the point at which
-- the value should be evaluated and the time of the current time-step.
-- Note that we use math-functions from luas standard math-library.
-- (here the . operator is used, since math is not an object but a library)
function inletVel2d(x, y, t)
	local Umax = 1.5
	return (1.0-y*y) * Umax, 0.0
--	return 1.0, -1.0
--	return 0.0, 0.0
end

function WallVel2d(x, y, t)
	return 0.0, 0.0
end

function outletVel2d(x, y, t)
	local Umax = 1.5
	return (1.0-y*y) * Umax
end

-- setup Outlet
OutletDisc = DirichletBoundary()
OutletDisc:add(0.0, "p", "Outlet")
OutletDisc:add(0.0, "v", "Outlet")

-- setup Inlet
InletDisc = NavierStokesInflow("u,v,p", "Inner")
InletDisc:add("inletVel"..dim.."d", "Inlet")

--setup Walles
WallDisc = NavierStokesWall("u,v,p")
WallDisc:add("UpperWall,LowerWall")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(InletDisc)
domainDisc:add(WallDisc)
domainDisc:add(OutletDisc)

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

-- In ug4 we use Operator interfaces. An operator is simply a mapping from in
-- space into the other. A frequently used implementation of such a mapping is
-- the usage of discretizations, that map a solution (grid function) onto some
-- right-hand side. We can create an operator that uses the recently created
-- domainDisc.
op = AssembledOperator(domainDisc)


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
domainDisc:adjust_solution(u)

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
--gmg:set_debug(dbgWriter)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(ConvCheck(100, 1e-16, 1e-8, true))
linSolver:set_compute_fresh_defect_when_finished(true)

-- choose a solver
solver = linSolver

-- Next we need a convergence check, that computes the defect within each
-- newton step and stops the iteration when a specified creterion is fullfilled.
-- For our purpose is the ConvCheck is sufficient. Please note,
-- that this class derives from a general IConvergenceCheck-Interface and
-- also more specialized or self-coded convergence checks could be used.
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(100)
newtonConvCheck:set_minimum_defect(1e-16)
newtonConvCheck:set_reduction(1e-6)
newtonConvCheck:set_verbose(true)

-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(5)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- Finally we set the non-linear operator created above and initiallize the
-- Newton solver for this operator.
newtonSolver:init(op)

-- In a first step we have to prepare the newton solver for the solution u. This
-- sets e.g. dirichlet values in u.
if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed."); exit();
end


SaveVectorForConnectionViewer(u, "StartSolution.vec")

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
vtkWriter:select_nodal("p", "pressure")
vtkWriter:print("Solution", u)