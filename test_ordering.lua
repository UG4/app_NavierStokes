-------------------------------------------------------------------------------------------------------
--
--  Driven cavity problem
--
--  Author: Christian Wehner
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.
elemType = util.GetParam("-elem", "quad")

InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 2 then
	if elemType == "tri" then 
		gridName = util.GetParam("-grid", "grids/dc_tri.ugx")
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	else
		gridName = util.GetParam("-grid", "grids/dc_quads.ugx")
		gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
	end
else print("Chosen Dimension " .. dim .. "not supported. Exiting."); exit(); end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",0)

print(" Chosen Parameters:")
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
requiredSubsets = {"Inner", "Top", "Bottom", "Right", "Left"}
requiredSubsets = {}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
--if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
approxSpace:add_fct("p", "piecewise-constant") 

-- finally we print some statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--OrderLex(approxSpace, "lr");
--OrderCuthillMcKee(approxSpace, true);

u = GridFunction(approxSpace)
u:set(0)

function StartValue_u(x,y,t) return 1 end
function StartValue_v(x,y,t) return 2 end
function StartValue_p(x,y,t) return 3 end

Interpolate("StartValue_u", u, "u")
--Interpolate("StartValue_v", u, "v")
--Interpolate("StartValue_p", u, "p")

SaveVectorForConnectionViewer(u, "Test.vec")
