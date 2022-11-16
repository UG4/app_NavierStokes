--------------------------------------------------------------------------------
--
--   Lua - Script to compute a channel problem
--
--	Navier-Stokes test problem on the unit square using dirichlet conditions
--  c.f. 
--		Nigon, P., Une nouvelle classe de methodes multigrilles pour 
--					les problemes mixtes, E.C.L. 84-19. Lyon 1984
--		Wittum, G., Multi-Grid Methods for Stokes and Navier-Stokes Equations,
--					Numer. Math. 54, 543-563, 1989
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("transforming.lua")

local dim 		= util.GetParamNumber("-dim", 2, "world dimension")
local numRefs 	= util.GetParamNumber("-numRefs", 3, "number of grid refinements")
local numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")

local order 		= util.GetParamNumber("-order", 2, "order pressure and velocity space")
local vorder 		= util.GetParamNumber("-vorder", order, "order velocity space")
local porder 		= util.GetParamNumber("-porder", order-1, "order pressure space")
local type     	= util.GetParam("-type", "fe", "Type of discretization")
local bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
local bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
local bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
local bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
local upwind      = util.GetParam("-upwind", "lps", "Upwind type")
local stab        = util.GetParam("-stab", "fields", "Stabilization type")
local diffLength  = util.GetParam("-difflength", "raw", "Diffusion length type")

-- Channel sizes as in grid file: [0,20] x [-1,1]
if 	dim == 2 then gridName = util.GetParam("-grid", "grids/unit_square_01_quads_2x2_pressure_node.ugx")
else print("Choosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

-- Lets write some info about the choosen parameter
print(" Selected Parameters:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. type)
--print("    only stokes      = " .. tostring(bStokes))
--print("    no laplace       = " .. tostring(bNoLaplace))
--print("    exact jacobian   = " .. tostring(bExactJac))
--print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

-- Init UG
InitUG(dim, AlgebraType("CPU", 1));

-- Create the domain and load a grid
local neededSubsets = {"Inner", "Boundary", "PressureNode"}
local dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- All subset are ok. So we can create the Approximation Space
local approxSpace = ApproximationSpace(dom)

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
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--OrderLex(approxSpace, "lr")
--OrderCuthillMcKee(approxSpace, false)  -- order vertices (crucial for)
--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

function Source2d(x, y, t) 
	return 36 * math.sin(3*(x+y)), 0				
end

-- create NavierStokes disc
NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_source("Source2d");
NavierStokesDisc:set_kinematic_viscosity(1);

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
function exactSolU2d(x, y, t) return math.sin(3*(x+y)) 				end
function exactSolV2d(x, y, t) return -math.sin(3*(x+y)) 			end
function exactSolP2d(x, y, t) return -6*math.cos(3*(x+y))			end
function exactSolVel2d(x, y, t)
--	return 0, 0
	return exactSolU2d(x,y,t), exactSolV2d(x,y,t)
end


FixPressureDisc = DirichletBoundary()
FixPressureDisc:add(0, "p", "PressureNode")--"Boundary, PressureNode"

BndDisc = NavierStokesInflow(NavierStokesDisc)
BndDisc:add("exactSolVel"..dim.."d", "Boundary, PressureNode")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(BndDisc)
--domainDisc:add(FixPressureDisc)

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------
u = GridFunction(approxSpace)
--w = GridFunction(approxSpace)

LUSolver = LU()
LUSolver:set_minimum_for_sparse(100000)

cmpGS = ComponentGaussSeidel(0.8, {"p"}, {0,1,2,1,0}, {1,1,1,1})
vankaGS = ElementGaussSeidel(0.9, "vertex")
local vankaSmoother = ComponentGaussSeidel(1.0, {"p"})
vankaSmoother:set_alpha(1.0) -- 1.0
vankaSmoother:set_beta(1.0)
vankaSmoother:set_weights(true)

transfer = StdTransfer()
--transfer:set_restriction_damping(0.5)

-- Linear Solver
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(2)
gmg:set_base_solver(LUSolver)

--gmg:set_smoother(ElementGaussSeidel(0.9, "vertex"))
gmg:set_smoother(vankaSmoother)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_cycle_type(1)
gmg:set_prolongation(transfer)

gmg:set_restriction(transfer)
--gmg:add_prolongation_post_process(AverageComponent("p"))
--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
gmg:set_rap(false)


cmpConvCheck = CompositeConvCheck(approxSpace, 10, 1e-8, 1e-8)
cmpConvCheck:set_component_check("u", 1e-8, 1e-8)
cmpConvCheck:set_component_check("v", 1e-8, 1e-8)
cmpConvCheck:set_component_check("p", 1e-8, 1e-8)

linSolver = LinearSolver()
--linSolver:set_preconditioner(ElementGaussSeidel(0.9, "vertex"))
linSolver:set_preconditioner(vankaSmoother)
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(ConvCheck(100, 1e-10, 1e-8, true))
linSolver:set_compute_fresh_defect_when_finished(true)
--linSolver:set_debug(GridFunctionDebugWriter(approxSpace))


-- Non-Linear Solver
newtonSolver = NewtonSolver(AssembledOperator(domainDisc))
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(1, 1e-8, 1e-6, true))
--newtonSolver:set_line_search(StandardLineSearch(5, 1, 0.5, true))
--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

vtkWriter = VTKOutput()
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "pressure")


local dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)
  
local dbgConvCheck = ConvCheck(20, 1e-10, 1e-8)
local dbgSmooth=u:clone()

local dbgSolver = LinearSolver()
dbgSolver:set_preconditioner(vankaSmoother) -- cgs, gmg
dbgSolver:set_convergence_check(dbgConvCheck)

local dbgIter= DebugIterator()
dbgIter:set_preconditioner(gmg)  -- gmg is the 'real' preconditioner
dbgIter:set_solver(dbgSolver)
dbgIter:set_solution(dbgSmooth)
dbgIter:set_random_bounds(-5e-6, 5e-6)
dbgIter:set_debug(dbgWriter)  -- print t_0 anf t_N



-- Interpolate Start Iterate
print("Interpolate")
Interpolate("exactSolU"..dim.."d", u, "u")
Interpolate("exactSolV"..dim.."d", u, "v")
Interpolate("exactSolP"..dim.."d", u, "p")
u:set(0.0)
--u:set_random(0,1)
vtkWriter:print("NigonStart", u)


A = MatrixOperator()
b = GridFunction(approxSpace) 
domainDisc:assemble_linear(A, b)
domainDisc:adjust_solution(u)
--[[
linSolver:set_debug(GridFunctionDebugWriter(approxSpace))

solverConvCheck = CompositeConvCheck(approxSpace, 100, 1e-12, 1e-20)
solverConvCheck:set_component_check({"u", "v", "p"}, 1e-12, 1e-20)

linSolver:set_convergence_check(solverConvCheck)
linSolver:init(A, u)
linSolver:apply(u,b)
vtkWriter:print("Nigon", u)
--exit()
--]]

-- Apply the newton solver. A newton itertation is performed to find the solution.
--if newtonSolver:apply(u) == false then
--	 print ("Newton solver apply failed."); exit();
--end





-- TODO: Adjust path
-- local loadtrafosmoothers = assert(loadfile("transforming.lua"))
--loadtrafosmoothers()







local trafoSmoother=util.transforming.CreateTransformingSmoother(approxSpace, VelCmp)
----local trafoSmoother=CreateSchurComplementSmoother(approxSpace, domainDisc, 2.0)
--local trafoSmoother=CreateBlockTrafoSmoother0(approxSpace, domainDisc, 3.0)
--local trafoSmoother=CreateBlockTrafoSmoother1(approxSpace, domainDisc,0.1)

--trafoSmoother:set_debug(dbgWriter)

--smoother = trafoSmoother
smoother = vankaGS
smoother = vankaSmoother
-- smoother = gmg

print("Smoother")
linSolver = LinearSolver()
linSolver:set_preconditioner(smoother)
stdConvCheck = ConvCheck(7, 1e-10, 1e-8, true)
linSolver:set_convergence_check(cmpConvCheck)
linSolver:set_compute_fresh_defect_when_finished(true)
linSolver:set_debug(dbgWriter)



--- run problem with plain smoothing

newtonSolver:set_linear_solver(linSolver)
linSolver:init(A, u)
u:set_random(-1.0, 1.0)
b:set(0)
linSolver:apply(u,b)
vtkWriter:print("NigonSolutionTSmooth", u)    -- print solution
--SaveMatrixForConnectionViewer(u, A, "Jacobian.mat") 
--quit()



-- run problem using MG solver
print("GMG")
gmg:set_smoother(smoother)
print (gmg:config_string())
linSolver:set_preconditioner(gmg)
linSolver:init(A, u)
u:set_random(-1.0, 1.0)
b:set(0)
linSolver:apply(u,b)
vtkWriter:print("NigonSolutionGMG", u)       -- print solution
