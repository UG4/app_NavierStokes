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

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
numRefs 	= util.GetParamNumber("-numRefs", 0, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")

order 		= util.GetParamNumber("-order", 1, "order pressure and velocity space")
vorder 		= util.GetParamNumber("-vorder", order, "order velocity space")
porder 		= util.GetParamNumber("-porder", order-1, "order pressure space")

type     	= util.GetParam("-type", "fv1", "Type of discretization")
bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type")
stab        = util.GetParam("-stab", "fields", "Stabilization type")
diffLength  = util.GetParam("-difflength", "raw", "Diffusion length type")

Umax 		= util.GetParamNumber("-umax", 1.5)
Viscosity   = util.GetParamNumber("-visco", 0.01)

-- Channel sizes as in grid file: [0,20] x [-1,1]
if 	dim == 2 then gridName = util.GetParam("-grid", "grids/channel20Q18x10.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- choose some solver:
sol = util.GetParam("-sol", "gmg")
numPreSmooth = util.GetParamNumber("-numPreSmooth", 2)
numPostSmooth = util.GetParamNumber("-numPostSmooth", 2)
if util.HasParamOption("-numSmooth") then
numPreSmooth = util.GetParamNumber("-numSmooth", 2)
numPostSmooth = numPreSmooth
end
smooth = util.GetParam("-smooth", "ilut", "Type of smoother [jac |Êilu | ilut |Êegs |Êgs |Êsgs]")
groupType = util.GetParam("-group", "element", "Type of grouping for ElementGaussSeidel (egs)")
baseLev = util.GetParamNumber("-baseLev", 0)
cycleType =  util.GetParam("-cycle", "V", "gmg-cycle type [V | W | F]")
bRAP = util.HasParamOption("-rap", "use rap product as level matrices")

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. type)
print("    only stokes      = " .. tostring(bStokes))
print("    no laplace       = " .. tostring(bNoLaplace))
print("    exact jacobian   = " .. tostring(bExactJac))
print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)
print("    Umax             = " .. Umax)
print("    Viscosity        = " .. Viscosity)
print("")
print(" solution process related parameters chosen:")
print("    linear solver   = " .. sol)
print("    gmg cycle       = " .. cycleType)
print("    gmg base level  = " .. baseLev)
print("    gmg smoother    = " .. smooth)
print("    # pre-smooth    = " .. numPreSmooth)
print("    # post-smooth   = " .. numPostSmooth)
print("    rap             = " .. tostring(bRAP))

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
if type == "fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
elseif type == "fv" then
	approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
	approxSpace:add_fct("p", "Lagrange", porder) 
elseif type == "fe" then
	if porder==0 then
		approxSpace:add_fct(VelCmp, "Crouzeix-Raviart",1)
		approxSpace:add_fct("p", "piecewise-constant") 
	else
		approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
elseif type=="fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else print("Disc Type '"..type.."' not supported."); exit(); end

-- finally we print some statistic on the distributed dofs
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

-- create NavierStokes disc
NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(Viscosity);

--upwind if available
if type == "fv1" or type == "fvcr" then
	NavierStokesDisc:set_upwind(upwind)
	NavierStokesDisc:set_peclet_blend(bPecletBlend)
end

-- fv1 must be stablilized
if type == "fv1" then
	NavierStokesDisc:set_stabilization(stab, diffLength)
end

-- fe must be stabilized for (Pk, Pk) space
if type == "fe" and porder == vorder then
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

baseLU = LU()
baseLU:set_minimum_for_sparse(10000)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_base_level(baseLev)
gmg:set_base_solver(baseLU)
gmg:set_gathered_base_solver_if_ambiguous(true)
if 	    smooth == "ilu"  then smoother = ILU();
elseif 	smooth == "ilut" then smoother = ILUT(1e-6);
elseif 	smooth == "egs"  then smoother = ElementGaussSeidel(groupType);
elseif 	smooth == "cgs"  then smoother = ComponentGaussSeidel(0.8, "p", {0,1,2,1,0}, {1,1,1,1})
elseif 	smooth == "jac"   then smoother = Jacobi(0.66);
elseif 	smooth == "gs"   then smoother = GaussSeidel();
elseif 	smooth == "sgs"  then smoother = SymmetricGaussSeidel();
else print("Smoother not set, use -smooth option [ilu, ilut, jac, egs, cgs, gs, sgs]"); exit(); end
gmg:set_smoother(smoother)
gmg:set_cycle_type(cycleType)
gmg:set_num_presmooth(numPreSmooth)
gmg:set_num_postsmooth(numPostSmooth)
gmg:set_rap(bRAP)

-- create Linear Solver
SmoothSolver = LinearSolver()
SmoothSolver:set_preconditioner(smoother)

-- create Linear Solver
GMGSolver = LinearSolver()
GMGSolver:set_preconditioner(gmg)

-- create BiCGStab Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)

-- select some of the created solver
if 		sol == "gmg" then linSolver = GMGSolver;
elseif 	sol == "bicgstab" then linSolver = BiCGStabSolver;
elseif 	sol == "smooth" then linSolver = SmoothSolver;
else print("Linear solver not set, use -sol option [gmg, bicgstab, smooth]"); exit(); end
linSolver:set_convergence_check(ConvCheck(60, 1e-10, 1e-8, true))

-- Non-Linear Solver
newtonSolver = NewtonSolver(AssembledOperator(domainDisc))
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(40, 1e-8, 1e-6, true))
--newtonSolver:set_line_search(StandardLineSearch(5, 1, 0.5, true))
--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

-- Interpolate Start Iterate
u = GridFunction(approxSpace)
u:set(0.0)

if util.HasParamOption("-startWithExact", "StartWithExactSol") == true then 
	Interpolate("exactSolU"..dim.."d", u, "u")
	Interpolate("exactSolV"..dim.."d", u, "v")
	Interpolate("exactSolP"..dim.."d", u, "p")
end

-- Apply the newton solver. A newton itertation is performed to find the solution.
if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

-- Output of solution
vtkWriter = VTKOutput()
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "pressure")
vtkWriter:print("Channel", u)