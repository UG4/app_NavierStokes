
-------------------------------------------------------------------------
--
--  Lua - Script to study the grid convergence of the 2d fixed cylinder
--  benchmark on a boundary fitted mesh. It applies the standard Navier
--  Stokes discretisation within the ug4 framework which was studied by
--  S. Naegele in her PhD thesis.
--  Experimental and computational reference data is provided by:
--
--  M. Schäfer, S. Turek: 'Benchmark computations of laminar flow around
--   'a cylinder'
--  in: Flow Simulation with High-Performance Computers II,
--  in: Notes on Numerical Fluid Mechanics, vol. 2, (1996), pp. 547–566.
--
--  and more recent reference values can be found in:
--
--  V. John, G. Matthies: 'Higher-order finite element discretizations in
--   'a benchmark problem for incompressible flows'
--  Int. J. Numer. Methods Fluids, vol. 37, (2001), pp. 885–903.
--
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
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
--  Setup Parameter
--------------------------------------------------------------------------------

local Viscosity	= 1e-3
local Um = 0.3
if dim == 3 then Um = 0.45 end
local H = 0.41
local L = 0.1
local Umean2 = math.pow(2/3*Um, 2)

lev = 6
base = 1

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/conv_rates_static.lua")

-- constants
dim		   = util.GetParamNumber("-dim", 2, "Dimension of the problem")
numPreRefs = util.GetParamNumber("-numPreRefs", base, "number of refinements")
numRefs    = util.GetParamNumber("-numRefs", lev, "number of refinements")

-- choose grid
if dim == 2 then
	gridName = util.GetParam("-grid", "grids/cylinder_tri_regular.ugx")
else
	print("Dimension "..dim.." not supported."); exit();
end


--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

function inletVel2d(x, y, t)
    return 4 * Um * y * (H-y) / (H*H), 0.0
end

function ZeroSolution(x, y, t)
	return 0.0
end

--------------------------------------------------------------------------------
--  Setup Domain
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", dim+1));
	
	local dom = Domain()
	LoadDomain(dom, gridName)
	
	
-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
	if     dim == 2 then
		ProjectVerticesToSphere(dom, {0.2, 0.2}, 0.05, 0.001)
        falloffProjector = SphereProjector(MakeVec(0.2, 0.2, 0), 0.1, 0.15)
	elseif dim == 3 then
        falloffProjector = CylinderProjector(MakeVec(0.5, 0.2, 0.0), MakeVec(0, 0, 1), 0.04, 0.1)
	end

	local projHandler = ProjectionHandler(dom:subset_handler())
    dom:set_refinement_projector(projHandler)

    projHandler:set_projector("Inner", falloffProjector)
    projHandler:set_projector("CylinderWall", falloffProjector)
    if dim == 3 then
        projHandler:set_projector("BackWall", falloffProjector)
        projHandler:set_projector("FrontWall", falloffProjector)
    end

-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
    local refiner =  GlobalDomainRefiner(dom)

	write("Pre-Refining("..numPreRefs.."): ")
	for i=1,numPreRefs do write(i .. " ");	refiner:refine(); end
	write("done. Distributing...")
    if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end
	write("done.\n")
	
	--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
	
	return dom
end


--------------------------------------------------------------------------------
--  Setup FV Element Discretization
--------------------------------------------------------------------------------

function CreateApproxSpace(dom, discType, p)
	local approxSpace = ApproximationSpace(dom)
	if dim >= 1 then approxSpace:add_fct("u", "Lagrange", 1) end
	if dim >= 2 then approxSpace:add_fct("v", "Lagrange", 1) end

	approxSpace:add_fct("p", "Lagrange", 1)
	return approxSpace
end


function CreateDomainDisc(approxSpace, discType, p)
	
	fctUsed = "u"
	fctUsed = fctUsed .. ", v" 
	fctUsed = fctUsed .. ", p"

    local stokes = true

    local NavierStokesDisc = NavierStokes(fctUsed, "Inner")
    stab = NavierStokesFIELDSStabilization()

    if not(stokes) then
        upwind = NavierStokesFullUpwind();
        stab:set_upwind(upwind)
    end
    stab:set_diffusion_length("RAW")
    NavierStokesDisc:set_stabilization(stab)
    if not(stokes) then
        NavierStokesDisc:set_upwind(upwind)
        NavierStokesDisc:set_pac_upwind(true)
    end

    NavierStokesDisc:set_peclet_blend(true)
    NavierStokesDisc:set_exact_jacobian(false)
	NavierStokesDisc:set_stokes(stokes)
	NavierStokesDisc:set_laplace(false)
	NavierStokesDisc:set_kinematic_viscosity(Viscosity);

--setup Inlet
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
--set pressure value in node
	PressureBndCond = DirichletBoundary()
	PressureBndCond:add("ZeroSolution", "p", "PressureBndCond")
 
--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	if dim == 2 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,PressureBndCond")
	elseif dim == 3 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall")	
	end
	
-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(PressureBndCond)
	domainDisc:add(InletDisc)
	domainDisc:add(WallDisc)
	
	return domainDisc
	
end

--------------------------------------------------------------------------------
--  Setup Solver
--------------------------------------------------------------------------------

function CreateSolver(approxSpace)

	dbgWriter = GridFunctionDebugWriter(approxSpace)
 
	local gmg = GeometricMultiGrid(approxSpace)
	gmg:set_base_solver(SuperLU())
	gmg:set_smoother(ILU())
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)

	local solver = SuperLU()
	
	local newtonConvCheck = ConvCheck()
	newtonConvCheck:set_maximum_steps(2000)
	newtonConvCheck:set_minimum_defect(1e-10)
	newtonConvCheck:set_reduction(1e-8)
	newtonConvCheck:set_verbose(true)
	
-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
	local newtonLineSearch = StandardLineSearch()
	newtonLineSearch:set_maximum_steps(10)
	newtonLineSearch:set_lambda_start(1.0)
	newtonLineSearch:set_reduce_factor(0.5)
	newtonLineSearch:set_accept_best(true)
	
-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(newtonConvCheck)
	newtonSolver:set_line_search(newtonLineSearch)
	if lev == base then
		newtonSolver:set_debug(dbgWriter)
	end
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)

    util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
    AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Compute Conv Rates
--------------------------------------------------------------------------------

options = {	

	size = 				{12.5, 6.75}, -- the size of canvas (i.e. plot)
	sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
	font = 				"Arial",
	fontsize =			12,
	
	logscale = 			true,
	grid = 				"lc rgb 'grey70' lt 0 lw 1", 
	linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
						linewidth = 3, pointsize = 1.3},
	border = 			" back lc rgb 'grey40' lw 2",
	decimalsign = 		",",
	key =	 			"on box lc rgb 'grey40' right top Left reverse spacing 2 width 1.1 samplen 2 height 0.5",
	tics =	 			{x = "nomirror out scale 0.75 format '%g' font ',8'",
						 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
	mtics =				5,
	slope = 			{dy = 3, quantum = 0.25, at = "last"},
	padrange = 			{ x = {0.8, 4}, y = {0.5, 1.5}},
	
	--multiplot = 		true, 
	--multiplot = 		{rows = 2, conjoined = true}
}

if util.HasParamOption("-replot") then
	util.rates.static.replot(options); exit()
end

			  
util.rates.static.compute(
{
	--ExactSol = {["u"] = "ExactSolution_u", ["v"] = "ExactSolution_v", ["p"] = "ExactSolution_p"},
	--ExactGrad =  {["c"] = "ExactGrad"..dim.."d"},
	
	CreateDomain = CreateDomain,
	CreateApproxSpace = CreateApproxSpace,
	CreateDomainDisc = CreateDomainDisc,
	CreateSolver = CreateSolver,
	
    ComputeSolution = ComputeNonLinearSolution,

	exact = false,
	maxlevel = true,
	prevlevel = true,
	interpol = false,
	
	DiscTypes = 
	{
	--{type = "fe", pmin = 1, pmax = 4, lmin = 0},
	{type = "fv1", pmin = 1, pmax = 1, lmin = numPreRefs},
	--{type = "fv1", pmin = 1, pmax = 1},
	--{type = "fvcr", pmin = 1, pmax = 1}
	},
	
	PrepareInitialGuess = function (u, lev, minLev, maxLev, domainDisc, solver)
			u[lev]:set(0.0)
	end,
		
	gpOptions = options,
	
	MaxLevelPadding = function(p)
 		return math.floor((p)/2)
	 end,

})
