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
 gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
 --gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
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

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Create the domain and load a grid
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we add the components of the velocity as Lagrange Ansatz functions of first order
--if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart", 1) end
--if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart", 1) end
--if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart", 1) end

-- we add the pressure as Lagrange Ansatz function of first order
approxSpace:add_fct("p", "Crouzeix-Raviart", 1) -- TODO: must be Constant Approx Space

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

u = GridFunction(approxSpace)

function StartValue_u(x,y,t) return x end
function StartValue_v(x,y,t) return y end
function StartValue_p(x,y,t) return x*y end

--Interpolate("StartValue_u", u, "u")
--Interpolate("StartValue_v", u, "v")
Interpolate("StartValue_p", u, "p")

SaveVectorForConnectionViewer(u, "StartSol.vec")

out = VTKOutput()
out:clear_selection()
out:select_all(false)
--out:select_element("u", "u")
--out:select_element("v", "v")
out:select_element("p", "p")
out:print("StartValue", u)
