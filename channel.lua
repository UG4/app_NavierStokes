--------------------------------------------------------------------------------
--
--   Lua - Script to compute a channel problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the channel problem
--
--   Author: Raphael Prohl, Andreas Vogel
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim 		= util.GetParamNumber("-dim", 2)
numRefs 	= util.GetParamNumber("-numRefs", 0)
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)

order 		= util.GetParamNumber("-order", 1)
vorder 		= util.GetParamNumber("-vorder", order)
porder 		= util.GetParamNumber("-porder", order-1)

discType   	= util.GetParam("-type", "fv1")
bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")

Umax 		= util.GetParamNumber("-umax", 1.5)
Viscosity   = util.GetParamNumber("-visco", 0.01)

-- Channel sizes as in grid file: [0,20] x [-1,1]
if 	dim == 2 then gridName = util.GetParam("-grid", "grids/channel20Q18x10.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

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
print("    Umax             = " .. Umax)
print("    Viscosity        = " .. Viscosity)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Init UG
InitUG(dim, AlgebraType("CPU", 1));

-- Create the domain and load a grid
neededSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"}; FctCmp = {"u", "p"};
elseif  dim == 2 then VelCmp = {"u", "v"}; FctCmp = {"u", "v", "p"};
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"};
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- we add the velocity and pressure as Lagrange Ansatz function of first order
if discType == "fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
elseif discType == "fv" then
	approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
	approxSpace:add_fct("p", "Lagrange", porder) 
elseif discType == "fe" then
	if porder==0 then
		approxSpace:add_fct(VelCmp, "Crouzeix-Raviart",1)
		approxSpace:add_fct("p", "piecewise-constant") 
	else
		approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
elseif discType=="fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else print("Disc Type '"..discType.."' not supported."); exit(); end

-- finally we print some statistic on the distributed dofs
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

-- create NavierStokes disc
NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(Viscosity);

--upwind if available
if type == "fv1" or type == "fvcr" then
	--upwind = NavierStokesNoUpwind();
	--upwind = NavierStokesFullUpwind();
	--upwind = NavierStokesSkewedUpwind();
	upwind = NavierStokesLinearProfileSkewedUpwind();
	--upwind = NavierStokesRegularUpwind();
	--upwind = NavierStokesPositiveUpwind();
	NavierStokesDisc:set_conv_upwind(upwind)
	
	NavierStokesDisc:set_peclet_blend(bPecletBlend)
end

-- fv1 must be stablilized
if type == "fv1" then
	stab = NavierStokesFIELDSStabilization()
	--stab = NavierStokesFLOWStabilization()
	stab:set_upwind(upwind)

	--stab:set_diffusion_length("RAW")
	stab:set_diffusion_length("FIVEPOINT")
	--stab:set_diffusion_length("COR")

	NavierStokesDisc:set_stabilization(stab)
end

-- fe must be stabilized for (Pk, Pk) space
if discType == "fe" and porder == vorder then
	NavierStokesDisc:set_stabilization(3)
end


-- exact solution of the problem
function exactSolU2d(x, y, t) return (1.0-y*y) * Umax 				end
function exactSolV2d(x, y, t) return 0.0              				end
function exactSolP2d(x, y, t) return -2 * Umax * Viscosity * (x-20) end
function exactSolVel2d(x, y, t)
	return exactSolU2d(x,y,t), exactSolV2d(x,y,t)
end


-- setup Outlet
OutletDisc = DirichletBoundary()
OutletDisc:add(0.0, "p", "Outlet")
OutletDisc:add(0.0, "v", "Outlet")

-- setup Inlet
InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("exactSolVel"..dim.."d", "Inlet")

--setup Walles
WallDisc = NavierStokesWall(NavierStokesDisc)
WallDisc:add("UpperWall,LowerWall")

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

-- Linear Solver
gmg = GeometricMultiGrid(approxSpace)
gmg:set_base_solver(LU())
gmg:set_smoother(ILUT())
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(100, 1e-10, 1e-8, true))

-- choose a solver
linSolver = gmgSolver

-- Non-Linear Solver
newtonSolver = NewtonSolver(AssembledOperator(domainDisc))
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(10, 1e-12, 1e-6, true))
--newtonSolver:set_line_search(StandardLineSearch(5, 1, 0.5, true))
--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- Interpolate Start Iterate
u = GridFunction(approxSpace)
Interpolate("exactSolU"..dim.."d", u, "u")
Interpolate("exactSolV"..dim.."d", u, "v")
Interpolate("exactSolP"..dim.."d", u, "p")
--u:set(0.0)

-- Apply the newton solver. A newton itertation is performed to find the solution.
if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

-- Output of solution
vtkWriter = VTKOutput()
vtkWriter:select_all(false)
vtkWriter:select_nodal(VelCmp, "velocity")
vtkWriter:select_nodal("p", "pressure")
vtkWriter:print("Channel", u)