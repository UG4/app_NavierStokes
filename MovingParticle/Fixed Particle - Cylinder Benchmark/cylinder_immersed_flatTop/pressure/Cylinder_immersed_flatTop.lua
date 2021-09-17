
-------------------------------------------------------------------------
--
--  Lua - Script to compute the 2d fixed cylinder benchmark on an unfitted
--  mesh with immersed circular interface for the embedded cylinder. It
--  applies the 'MovingParticle' class as Navier Stokes discretisation
--  implemented to handle moving particles. For this application velocity
--  of the particle (= cylinder) is fixed to zero and therefore not part of
--  the system of unknowns.
--  The central output to be studied here is the profile of the pressure
--  along the immersed interface of the cylinder, see lua calls
--    movingParticle:print_pressure_nodal(u, lev)
--  and
--    movingParticle:movingParticle:print_deltaP(u, lev).
--  The adapted handling of the interface within the 'MovingParticle' class
--  avoid spurious oscillations of the pressure.
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
    gridName = util.GetParam("-grid", "grids/channel_41x220_tri_regular.ugx")
else print("Dimension "..dim.." not supported."); exit(); end


------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
-- Parameter Setup
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

-- particle radius and center
radius = 0.05
center = {0.2, 0.2}

-- viscosity of the fluid
visc = 1e-3

-- densities
densityPrt = 2.0
densityFluid = 1.0


StdFVAssembling = false
-- = true: uses non-projected, i.e. non-Flat-Top ansatz- and test-spaces
--          --> causes oscillations of the pressure along the cylinder interface
-- = false: removes the oscillations of the pressure along the cylinder interface

print_cutElemData = true
-- = true: prints computed coordinates of the corners of cut elements and local indices to a file

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

-- choose solver parameter
Solver = 1
newtonRed =  util.GetParamNumber("-newtonRed", 1e-12)


--if numPreRefs == 6 then
--numPreRefs = 5
--end

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

function CreateDomainParallel(neededSubsets)

 	local dom = Domain()
	LoadDomain(dom, gridName)
				
	refiner = GlobalDomainRefiner(dom)
	
	-- perform pre-refinement
	for i = 1, numPreRefs do
		print("refining...")
		refiner:refine()
	end

	partitionMap = PartitionMap()
	util.PartitionMapRegularGrid(dom, partitionMap, np_x, np_y, 1)

	if DistributeDomain(dom, partitionMap, true) == false then
		print("Distribution failed. Please check your partitionMap.")
		exit()
	end
	
	if util.CheckSubsets(dom, neededSubsets) == false then 
		print("Something wrong with required subsets. Aborting.");
		exit();
	end
	
	for i=numPreRefs+1,numRefs do
		refiner:refine()
	end

	delete(refiner)
	SaveDomain(dom, "distributed_domain_p" .. ProcRank() .. ".ugx")
	
	return dom
	
end


	neededSubsets = {"Inner", "Inlet", "Outlet", "Upper", "Lower", "PressureBndCond"}
	if NumProcs() == 1 then 
		dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)		
	else
		dom = CreateDomainParallel(neededSubsets)
	end
	

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

	NavierStokesDisc = NavierStokesFV1_cutElem(fctUsed, "Inner")

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
-- Set up fluid-particle problem
--------------------------------------------------------------------------
--------------------------------------------------------------------------

    domainDisc = DomainDiscretization(approxSpace)
    domainDisc:add(NavierStokesDisc)

--------------------------------------------------------------------------
-- Set up Particle functionality
--------------------------------------------------------------------------

-- add the cylinder to the discretisation
    interfaceProvider = ParticleProviderSphere()
    interfaceProvider:add_moving(radius, MakeVec(center[1], center[2]), densityPrt,
        MakeVec(0.0, 0.0), MakeVec(0.0, 0.0))

    cutElementHandler = CutElementHandler_FlatTop(dom:grid(), "u,v,p", interfaceProvider)

-- create the class handling the immersed particle
    movingParticle = MovingParticle(domainDisc, NavierStokesDisc, cutElementHandler, densityFluid, visc)

-- for the fixed cylinder, no gravity force is acting
    movingParticle:set_gravity(false, 0.0)

-- StdFVAssembling = false ==> removes the oscillations of the pressure along the cylinder interface
    movingParticle:set_StdFV_assembling(StdFVAssembling)

if SkipNearInterfaceElem == 1 then
    -- sets the threshold to '0.25 * meanLength*meanLength', with meanLenght = MeanElementDiameter(domain, lev)
    movingParticle:initialize_threshold(dom, numRefs, numRefs)

    -- it is also possible, to set a user defined threshold = 'thres' by calling
    --   movingParticle:set_threshold(lev, thres)
end

-- prints computed coordinates of all corners of cut elements and their local indices to a file
-- see './CutElementData.txt'
    movingParticle:set_print_cutElemData(print_cutElemData)


    interfaceProvider:print()


--------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------

--setup inlet and outlet
    InletDisc = NavierStokesInflowFV1_cutElem(NavierStokesDisc)
    InletDisc:add("inletVel2d", "Inlet, Outlet")

--setup Walles
    WallDisc = NavierStokesWall(NavierStokesDisc)
    WallDisc:add("Upper,Lower")

-- set coupling condition on immersed boundary
    InnerBndCond = movingParticle:get_BndCond()

--set pressure value in node
    PressureBndCond = DirichletBoundary()
    PressureBndCond:add("ZeroSolution", "p", "PressureBndCond")

-- add boundary conditions
    domainDisc:add(InnerBndCond)
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

-- we have to prepare the particle class. This sets the marker for the
-- cut elements and initializes the values
-- initialize 'movingIParticle' class with topLevel as argument TWICE:
--  ==> no loop for multigrid
    movingParticle:init(u, approxSpace, base, lev)

--  print the number of elements cut by the immersed cylinder line
    prtIndex = 0
    numCutElem = movingParticle:get_numCutElements(lev, prtIndex)
    print("\n #cut elements: " .. numCutElem)
    print("\n For detailed info on the cut element geometry see file ./'CutElementData.txt'\n")


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

-- write the pressure along the cylinder interface into file
    movingParticle:print_pressure_nodal(u, lev)

-- compute the pressure at the front and back of the cylinder; write it to a file
    movingParticle:print_deltaP(u, lev)

-- write start solution to ConnectionViewer
    SaveVectorForConnectionViewer(u, "Solution_cylinder_immersed_flatTop.vec")

-- write start solution to vtk
    vtkWriter = VTKOutput()
    vtkWriter:select_all(false)
    vtkWriter:select_nodal("u,v", "velocity")
    vtkWriter:select_nodal("u", "u")
    vtkWriter:select_nodal("v", "v")
    vtkWriter:select_nodal("p", "pressure")
    vtkWriter:print("Solution_cylinder_immersed_flatTop", u)

    print("done.")

