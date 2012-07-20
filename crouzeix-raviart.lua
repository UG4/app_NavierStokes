------------------------------------------------------------------
--
--   Lua - Script to test the Navier-Stokes implementation
--
--	This script sets up a problem for the Navier-Stokes discretization
--	on staggered grid using Crouzeix-Raviart type elements
--  (see cylinder.lua for collocated grid stabilization method)
--
--   Author: Christian Wehner
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
InitUG(dim, AlgebraType("CPU", 1));

-- Next, we decide which grid to use. This can again be passed as a command line
-- option or a default value is used.
if 	dim == 2 then
-- gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
gridName = util.GetParam("-grid", "grids/cylinder.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
neededSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall", "CylinderWall"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we add the components of the velocity as Lagrange Ansatz functions of first order
if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end

-- we add the pressure as Lagrange Ansatz function of first order
approxSpace:add_fct("p", "piecewise-constant") 

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

-- We set up the discretization. The basic idea is to first create the single
-- components of the discretization (e.g. element discretization, boundary
-- conditions, ...) and then glue everything together in one DomainDiscretization
-- object that does nothing else than grouping the single components. When it
-- comes to assembling this object runs all scheduled discretization component
-- one after the other.

-- Lets begin with the element discretization. This is the core part for the
-- Navier-Stokes equation and provides the local stiffness- and mass matrices,
-- that are then added to the global matrices. First we concatenate a string
-- that holds all function we use. This string must be passed to the navier
-- stokes elem disc to indicate which function is used for velocity and which
-- one for pressure. E.g. for dim==2 this gives: "u, v, p".
fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

-- Now we create the Navier-Stokes disc. Here we use a util function, that
-- creates the disc for the right world dimension according to the loaded domain
-- (this is handled internally, since the approxSpace knows the domain and the
-- domain knows the dimension). We tell the elem Disc to assemble on all elements
-- that are in the subset "Inner". If we would have several domains, where we
-- would like to do the same, this could be done by passing a list of subsets
-- separated by ',', (e.g. "Inner1, Inner2, Inner3").
elemDisc = NavierStokes(fctUsed, "Inner")

elemDisc:set_disc_scheme("staggered")

-- Now, we have to setup the stabilization, that is used for the Continuity Equation.
-- The stabilization is passed to the Navier-Stokes elem disc as an object.
-- Therefore, we created it now and will pass it to the disc. But first, we have
-- to create the Upwind scheme, that is used inside the stabilization. There are
-- several possibilities:

noUpwind = NavierStokesCRNoUpwind()
fullUpwind = NavierStokesCRFullUpwind()

-- ... and set the upwind
elemDisc:set_conv_upwind(fullUpwind)

elemDisc:set_peclet_blend(false)
elemDisc:set_exact_jacobian(false)
elemDisc:set_stokes(false)
elemDisc:set_laplace(false)

-- ... and finally we choose a value for the kinematic viscosity.
elemDisc:set_kinematic_viscosity(1.0e-3)


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

function u_inletVel2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return 4 * Um * y * (H-y) 
end

-- setup Outlet
dirichletBnd = DirichletBoundary()
dirichletBnd:add(0.0, "p", "Outlet")

-- setup Inlet
-- LuaInletDisc = NavierStokesInflow("u,v,p", "Inner")
-- LuaInletDisc:add("inletVel"..dim.."d", "Inlet")

-- setup Wall
-- WallDisc = NavierStokesWall("u,v,p")
-- WallDisc:add("UpperWall,LowerWall,CylinderWall")
dirichletBnd:add(0.0, "u", "UpperWall,LowerWall,CylinderWall")
dirichletBnd:add(0.0, "v", "UpperWall,LowerWall,CylinderWall")

uInletVel = LuaUserNumber("u_inletVel"..dim.."d")
dirichletBnd:add(uInletVel, "u", "Inlet")
dirichletBnd:add(0.0, "v", "Inlet")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
-- domainDisc:add(LuaInletDisc)
-- domainDisc:add(WallDisc)
domainDisc:add(dirichletBnd)

--------------------------------
--------------------------------
-- Operator
--------------------------------
--------------------------------

-- In ug4 we use Operator interfaces. An operator is simply a mapping from in
-- space into the other. A frequently used implementation of such a mapping is
-- the usage of discretizations, that map a solution (grid function) onto some
-- right-hand side. We can create an operator that uses the recently created
-- domainDisc.
op = AssembledOperator(domainDisc)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------

-- Now lets solve the problem. Create a vector of unknowns and a vector
-- which contains the right hand side. We will use the approximation space
-- to create the vectors. Make sure to create the vector for the same
-- dofs as set to linOp through linOp:set_dof_distribution.
-- Note that the vectors automatically have the right size.
u = GridFunction(approxSpace)

-- Init the vector representing the unknowns with 0 to obtain
-- reproducable results.
u:set(0)

function StartValue_u(x,y,t) return x end
function StartValue_v(x,y,t) return y end
function StartValue_p(x,y,t) return x end

Interpolate("StartValue_u", u, "u")
Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

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
linSolver = BiCGStab()
-- linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(StandardConvergenceCheck(1000000, 1e-7, 1e-12, true))

-- choose a solver
solver = linSolver
-- solver = LU()

-- Next we need a convergence check, that computes the defect within each
-- newton step and stops the iteration when a specified creterion is fullfilled.
-- For our purpose is the StandardConvergenceCheck is sufficient. Please note,
-- that this class derives from a general IConvergenceCheck-Interface and
-- also more specialized or self-coded convergence checks could be used.
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(5e-4)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)

-- Sometimes its helpful to write the defect and jacobian of the newton step
-- to debug the implementation. For that, we use the debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

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

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("CRSolution", u)

print("done.")
