
-------------------------------------------------------------------------
--
--  Lua - Script to compute the 2d fixed cylinder benchmark on a boundary
--  fitted mesh. It applies the standard Navier Stokes discretisation
--  within the ug4 framework which was studied by S. Naegele in her PhD
--  thesis.
--  The central output to be studied here is the difference of the pressure
--  between front and back of the immersed cylinder for comparison with
--  the unfitted approach applying the 'MovingParticle' class.
--
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
-----------------------------------------------------------------------------------


------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
-- Setup
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim         = util.GetParamNumber("-dim", 2, "Dimension of the problem")
numPreRefs  = util.GetParamNumber("-numPreRefs", 2, "number of refinements")
numRefs     = util.GetParamNumber("-numRefs", 2, "number of refinements")

InitUG(dim, AlgebraType("CPU", dim+1));

if dim == 2 then
    gridName = util.GetParam("-grid", "grids/cylinder_tri_regular.ugx")
else print("Dimension "..dim.." not supported."); exit(); end


------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
-- Parameter Setup
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
-- viscosity of the fluid
visc = 1e-3


-- choose solver parameter
Solver = 0
newtonRed =  util.GetParamNumber("-newtonRed", 1e-12)


if numPreRefs == 6 then
numPreRefs = 5
end

base = numPreRefs
lev = numRefs

-- parameter for parallel computations
np_x = 4
np_y = 1

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
print("	   visc             = " .. visc)
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

function CreateDomain()

    local dom = Domain()
    LoadDomain(dom, gridName)

-- NOTE: Projector creation in script-code is deprecated. Instead one should
--		 add projectors to individual subsets directly in ProMesh.
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

    for i=1,numPreRefs do write(i .. " ");	refiner:refine(); end
    if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
        exit();
    end
    for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end

    --SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)

    return dom
end

    dom = CreateDomain()

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

    Um = 0.3
    H = 0.41
    L = 0.1
    Umean2 = math.pow(2/3*Um, 2)

    function StartValInletVel2d(x, y, t)
        return 4 * Um * y * (H-y) / (H*H)
    end

    function inletVel2d(x, y, t)
        return 4 * Um * y * (H-y) / (H*H), 0.0
    end

    function ZeroSolution(x, y, t)
        return 0.0
    end

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Navier-Stokes Discretization Setup
--------------------------------------------------------------------------
--------------------------------------------------------------------------

	fctUsed = "u"
	fctUsed = fctUsed .. ", v" 
 	fctUsed = fctUsed .. ", p"

	NavierStokesDisc = NavierStokes(fctUsed, "Inner")

-- set upwind
    upwind = NavierStokesLinearProfileSkewedUpwind();

-- set stabilisation
    stab = NavierStokesFLOWStabilization()
    stab:set_upwind(upwind)

	--stab:set_diffusion_length("FIVEPOINT")
	--stab:set_diffusion_length("COR")
	stab:set_diffusion_length("RAW")

-- set up
	NavierStokesDisc:set_stabilization(stab)
    NavierStokesDisc:set_upwind(upwind)
    NavierStokesDisc:set_pac_upwind(true)
 	NavierStokesDisc:set_peclet_blend(true)
	NavierStokesDisc:set_exact_jacobian(false)
	NavierStokesDisc:set_stokes(false)
	NavierStokesDisc:set_laplace(false)
	NavierStokesDisc:set_kinematic_viscosity(visc);

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Set up problem
--------------------------------------------------------------------------
--------------------------------------------------------------------------

    domainDisc = DomainDiscretization(approxSpace)
    domainDisc:add(NavierStokesDisc)


--------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------

--setup inlet and outlet
    InletDisc = NavierStokesInflow(NavierStokesDisc)
    InletDisc:add("inletVel2d", "Inlet, Outlet")

--setup Walles
    WallDisc = NavierStokesWall(NavierStokesDisc)
    WallDisc:add("UpperWall,LowerWall,CylinderWall,PressureBndCond")

--set pressure value in node
    PressureBndCond = DirichletBoundary()
    PressureBndCond:add("ZeroSolution", "p", "PressureBndCond")

-- add boundary conditions
    domainDisc:add(InletDisc)
    domainDisc:add(WallDisc)
    domainDisc:add(PressureBndCond)

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Set up Solvers
--------------------------------------------------------------------------
--------------------------------------------------------------------------

    ilu = ILU()
	ilu:set_beta(0.99)

    sgs = SymmetricGaussSeidel()

    Base = LinearSolver()
	Base:set_convergence_check(ConvCheck(10000, 1e-13, 1e-13, true))
	Base:set_preconditioner(sgs)
 
    baseSolver = BiCGStab()
	baseSolver:set_preconditioner(ilu)
	baseSolver:set_convergence_check(ConvCheck(5000, 1e-10, 1e-10, false))
	baseSolver:set_compute_fresh_defect_when_finished(true)
  

    gmg = GeometricMultiGrid(approxSpace)
	gmg:set_base_level(numPreRefs)
	gmg:set_base_solver(SuperLU())
	gmg:set_smoother(ilu)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(12)
	gmg:set_num_postsmooth(12)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	linSolver = LinearSolver()	
	linSolver:set_preconditioner(gmg)
	linSolver:set_convergence_check(ConvCheck(10000, 1e-12, 1e-30, true))
	
	if do_transfer == true then
        transfer = ParticleTransfer(approxSpace, cutElementHandler)
		gmg:set_transfer(transfer)
		transfer:set_debug(dbgWriter)
		transfer:set_use_transposed(true)
    end

-- choose a solver
    if Solver == 1 then solver = SuperLU()
    else                solver = linSolver
    end

    newtonConvCheck = ConvCheck()
	newtonConvCheck:set_maximum_steps(2000)
	newtonConvCheck:set_minimum_defect(1e-12)
	newtonConvCheck:set_reduction(1e-12)
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
    u = GridFunction(approxSpace)
    u:set(0)

    Interpolate("StartValInletVel2d", u, "u");

    domainDisc:adjust_solution(u)

--------------------------------------------------------------------------
--------------------------------------------------------------------------
-- Init, prepare, apply newton Solver
--------------------------------------------------------------------------
--------------------------------------------------------------------------

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


-- compute the pressure at the front and back of the cylinder;
    PEval = GlobalGridFunctionNumberData(u, "p")
    pFront = PEval:evaluate_global({0.15, 0.2})
    pBack = PEval:evaluate_global( {0.25, 0.2} )
    Delta_P = pFront-pBack

-- write it to a file
    filename = "./delta_p_level"..numRefs..".txt"
    file = io.open(filename, "w")
    io.output(file)
    io.write(pFront.." \t "..pBack.." \t "..Delta_P.." # pFront, pBack, deltaP\n")
    io.close(file)

-- write start solution to ConnectionViewer
    SaveVectorForConnectionViewer(u, "Solution_cylinder_fitted.vec")

-- write start solution to vtk
    vtkWriter = VTKOutput()
    vtkWriter:select_all(false)
    vtkWriter:select_nodal("u,v", "velocity")
    vtkWriter:select_nodal("u", "u")
    vtkWriter:select_nodal("v", "v")
    vtkWriter:select_nodal("p", "pressure")
    vtkWriter:print("Solution_cylinder_fitted", u)

    print("done.")

