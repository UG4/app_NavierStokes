-------------------------------------------------------------------------------------------------------
--
--  Lua - Script to test the Navier-Stokes implementation
--
--  Author: Christian Wehner
--
--	A theoretical example to test the Navier-Stokes discretization.
--	The boundary conditions are inflow boundary conditions (Dirichlet conditions for the velocity)
--  on the whole boundary.
--  The analytical solution can be constructed e.g. via maple:
--
--  # construct u and v solution so that velocity is divergence-free
--  g:=x*x*x+x*y+y*y
--  u:=diff(g,y)
--  v:=-diff(g,x) 
--  # chose p
--  p:=x
--  # rhs is chosen so that Navier-Stokes system is fulfilled
--  rhsu:=factor(simplify(-1/R*diff(diff(u,x),x)-1/R*diff(diff(u,y),y)+u*diff(u,x)+v*diff(u,y)+diff(p,x)))
--  rhsv:=factor(simplify(-1/R*diff(diff(v,x),x)-1/R*diff(diff(v,y),y)+u*diff(v,x)+v*diff(v,y)+diff(p,y)))
--
------------------------------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- chose "fvcr" or "stabil"
discType = util.GetParam("-type", "fvcr")
order = util.GetParamNumber("-order", 1)
vorder = util.GetParamNumber("-vorder", order)
porder = util.GetParamNumber("-porder", order-1)
R=1

InitUG(dim, AlgebraType("CPU", 1))

if 	dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_1x1.ugx")
else
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_tets.ugx")
end
if dim~=2 and dim~=3 then
   print("Chosen Dimension " .. dim .. " not supported. Exiting.") exit() 
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs",2)

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. discType)
print("    grid       	= " .. gridName)
print("    v ansatz order = " ..vorder)
print("    p ansatz order = " ..porder)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Lets define a list of all subsets that we need
requiredSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

-- All subset are ok. So we can create the Approximation Space
approxSpace = ApproximationSpace(dom)

-- we destinguish the components of the velocity 
if 		dim == 1 then VelCmp = {"u"} FctCmp = {"u", "p"}
elseif  dim == 2 then VelCmp = {"u", "v"} FctCmp = {"u", "v", "p"}
elseif  dim == 3 then VelCmp = {"u", "v", "w"} FctCmp = {"u", "v", "w", "p"}
else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end

if discType=="fvcr" then
	if dim >= 1 then approxSpace:add_fct("u", "Crouzeix-Raviart") end
	if dim >= 2 then approxSpace:add_fct("v", "Crouzeix-Raviart") end
	if dim >= 3 then approxSpace:add_fct("w", "Crouzeix-Raviart") end
	approxSpace:add_fct("p", "piecewise-constant") 
end
if discType=="fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
end
if discType=="fv" or discType=="fe" then
	if porder==0 and vorder==1 then
		approxSpace:add_fct(VelCmp, "Crouzeix-Raviart",1)
	else
		approxSpace:add_fct(VelCmp, "Lagrange",vorder)
	end
	if porder==0 then
		approxSpace:add_fct("p", "piecewise-constant") 
    else 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
end

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

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", discType)

if discType=="fvcr" then
	-- set upwind
	noUpwind = NavierStokesNoUpwind()
	fullUpwind = NavierStokesFullUpwind()
	NavierStokesDisc:set_upwind(fullUpwind)
	NavierStokesDisc:set_defect_upwind(false)
	NavierStokesDisc:set_defect_upwind(true)
end	
if discType=="fv1" then
	upwind = NavierStokesNoUpwind()
	-- upwind = NavierStokesFullUpwind()
	-- upwind = NavierStokesSkewedUpwind()
	-- upwind = NavierStokesLinearProfileSkewedUpwind()
	--upwind = NavierStokesRegularUpwind()
	--upwind = NavierStokesPositiveUpwind()
	-- chose stabilization 
	stab = NavierStokesFIELDSStabilization()
	stab = NavierStokesFLOWStabilization()
	-- ... and set the upwind
	stab:set_upwind(upwind)
	--stab:set_diffusion_length("RAW")
	--stab:set_diffusion_length("FIVEPOINT")
	stab:set_diffusion_length("COR")
	
	-- set stabilization
	NavierStokesDisc:set_stabilization(stab)
	
	-- set upwind
	NavierStokesDisc:set_upwind(upwind)

end

NavierStokesDisc:set_peclet_blend(false)
NavierStokesDisc:set_exact_jacobian(false)
NavierStokesDisc:set_stokes(true)
NavierStokesDisc:set_laplace(true)
NavierStokesDisc:set_kinematic_viscosity(1.0/R)



----------------------------------
----------------------------------
-- Boundary conditions
----------------------------------
----------------------------------

function usol2d(x, y, t)
--	return 2*x^2*(1-x)^2*(y*(1-y)^2-y^2*(1-y))		
	return y
end

function vsol2d(x,y,t)
--	return -2*y^2*(1-y)^2*(x*(1-x)^2-x^2*(1-x))
	return -x
end

function psol2d(x,y,t)
--	return x^3+y^3+0.5
	return 0
--	return x+y
end

function usol3d(x,y,z,t)
--	return 2*x^2*y*z*(2*z-1)*(z-1)*(2*y-1)*(y-1)*(x-1)^2
	return y
end

function vsol3d(x,y,z,t)
--	return -x*y^2*z*(2*z-1)*(z-1)*(y-1)^2*(2*x-1)*(x-1)
	return -x
end

function wsol3d(x,y,z,t)
--	return -x*y*z^2*(z-1)^2*(2*y-1)*(y-1)*(2*x-1)*(x-1)
	return 0
end

function psol3d(x,y,z,t)
	return 0
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

function inletVel3d(x, y,z, t)
	return usol3d(x,y,z,t),vsol3d(x,y,z,t), wsol3d(x,y,z,t)
end

uSolution = LuaUserNumber("usol"..dim.."d")
vSolution = LuaUserNumber("vsol"..dim.."d")
if dim==3 then
	wSolution = LuaUserNumber("wsol"..dim.."d")
end
pSolution = LuaUserNumber("psol"..dim.."d")

InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Boundary")

----------------------------------
----------------------------------
-- Source
----------------------------------
----------------------------------

function source2d(x, y, t)
	return 
--	3*x*x+ -- dx p term
--	-4*(-1+2*y)*(-y+6*y*x-6*y*x^2+y^2-6*y^2*x+6*y^2*x^2+3*x^2-6*x^3+3*x^4)/R
	0
	,
--	3*y*y+ -- dy p term
--	4*(-1+2*x)*(3*y^2-6*y^3+3*y^4-x+6*y*x-6*y^2*x+x^2-6*y*x^2+6*y^2*x^2)/R
	0
end

function source3d(x, y,z, t)
	return 
--	3*x*x+ -- dx p term
-- 2*(-18*x*y*z^2-12*x^4*y*z^3+18*x^4*y*z^2+24*x^3*y*z^3-36*x^3*y*z^2-24*x^2*y*z^3+36*x^2*y^2*z-18*x*y^2*z-24*x^2*y^3*z^3+36*x^2*y^3*z^2+36*x^2*y^2*z^3-18*x^2*y*z+6*x*y*z+24*x^3*y*z-12*x^4*y*z+54*x*y^2*z^2-12*x^4*y^3*z+18*x^4*y^2*z+24*x^3*y^3*z-36*x^3*y^2*z-24*x^2*y^3*z+6*x^2*z^3-2*y^3*z+12*x*y^3*z+6*y^2*z^3+6*y^3*z^2-4*y^3*z^3-36*x*y^2*z^3-36*x*y^3*z^2+24*x*y^3*z^3-2*y*z^3+12*x*y*z^3-1440*x^4*y^5*z^5*R+1620*x^4*y^5*z^4*R+168*x^3*y^3*z^5*R-840*x^4*y^3*z^5*R+1620*x^4*y^4*z^5*R-54*x^2*y^2*z^2+36*x^5*y^2*z^2*R-20*x^4*y^2*z^2*R-20*x^3*y^3*z^2*R-20*x^3*y^2*z^3*R-160*x^4*y^6*z^6*R+940*x^4*y^3*z^4*R+336*x^7*y^3*z^5*R-376*x^7*y^3*z^4*R-1176*x^6*y^3*z^5*R+1316*x^6*y^3*z^4*R+1512*x^5*y^3*z^5*R-1692*x^5*y^4*z^3*R+940*x^4*y^4*z^3*R+576*x^7*y^5*z^5*R-648*x^7*y^5*z^4*R-648*x^7*y^4*z^5*R+728*x^7*y^4*z^4*R-2016*x^6*y^5*z^5*R+2268*x^6*y^5*z^4*R+2268*x^6*y^4*z^5*R-2548*x^6*y^4*z^4*R+2592*x^5*y^5*z^5*R-2916*x^5*y^5*z^4*R-2916*x^5*y^4*z^5*R+864*x^5*y^3*z^3*R-480*x^4*y^3*z^3*R-672*x^6*y^3*z^3*R+192*x^7*y^3*z^3*R-1820*x^4*y^4*z^4*R+336*x^7*y^5*z^3*R-376*x^7*y^4*z^3*R-1176*x^6*y^5*z^3*R+1316*x^6*y^4*z^3*R+1512*x^5*y^5*z^3*R-324*x^5*y^2*z^5*R+168*x^3*y^5*z^3*R-840*x^4*y^5*z^3*R-324*x^3*y^4*z^5*R-324*x^3*y^5*z^4*R-6*x^3*z+3*x^4*z-12*x^3*z^3+6*x^4*z^3+18*x^3*z^2-9*x^4*z^2+6*x^2*y^3-12*x^3*y^3+6*x^4*y^3+18*x^3*y^2-9*x^4*y^2-6*x^3*y+3*x^4*y+288*x^3*y^5*z^5*R+4*x^3*y^2*z^2*R+180*x^4*y^2*z^5*R-864*x^5*y^6*z^5*R+672*x^6*y^6*z^5*R-756*x^6*y^6*z^4*R+280*x^4*y^6*z^3*R+392*x^6*y^6*z^3*R-504*x^5*y^6*z^3*R+480*x^4*y^6*z^5*R-540*x^4*y^6*z^4*R+40*x^3*y^4*z^2*R+40*x^3*y^2*z^4*R-56*x^3*y^6*z^3*R+108*x^3*y^6*z^4*R-96*x^3*y^6*z^5*R+972*x^5*y^6*z^4*R+3276*x^5*y^4*z^4*R+108*x^5*y^2*z^6*R-60*x^4*y^2*z^6*R+140*x^6*y^2*z^3*R-40*x^7*y^2*z^3*R+252*x^6*y^2*z^5*R-72*x^7*y^2*z^5*R-280*x^6*y^2*z^4*R+80*x^7*y^2*z^4*R-324*x^5*y^5*z^2*R+252*x^6*y^5*z^2*R-72*x^7*y^5*z^2*R-280*x^6*y^4*z^2*R+80*x^7*y^4*z^2*R+140*x^6*y^3*z^2*R-40*x^7*y^3*z^2*R-36*x^3*y^5*z^2*R+972*x^5*y^4*z^6*R-864*x^5*y^5*z^6*R-756*x^6*y^4*z^6*R+672*x^6*y^5*z^6*R-504*x^5*y^3*z^6*R+392*x^6*y^3*z^6*R+288*x^5*y^6*z^6*R-224*x^6*y^6*z^6*R+12*x^3*y^6*z^2*R-28*x^6*y^2*z^2*R+32*x^3*y^6*z^6*R+180*x^4*y^5*z^2*R-60*x^4*y^6*z^2*R+108*x^5*y^6*z^2*R-84*x^6*y^6*z^2*R-84*x^6*y^2*z^6*R+8*x^7*y^2*z^2*R-192*x^7*y^6*z^5*R+216*x^7*y^6*z^4*R-112*x^7*y^6*z^3*R+216*x^7*y^4*z^6*R-192*x^7*y^5*z^6*R-112*x^7*y^3*z^6*R+64*x^7*y^6*z^6*R+24*x^7*y^6*z^2*R+24*x^7*y^2*z^6*R-1692*x^5*y^3*z^4*R-188*x^3*y^3*z^4*R-200*x^4*y^2*z^4*R+360*x^5*y^2*z^4*R+364*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+100*x^4*y^3*z^2*R-180*x^5*y^3*z^2*R+96*x^3*y^3*z^3*R+100*x^4*y^2*z^3*R-180*x^5*y^2*z^3*R-200*x^4*y^4*z^2*R+360*x^5*y^4*z^2*R-188*x^3*y^4*z^3*R+12*x^3*y^2*z^6*R-540*x^4*y^4*z^6*R+480*x^4*y^5*z^6*R+280*x^4*y^3*z^6*R+108*x^3*y^4*z^6*R-96*x^3*y^5*z^6*R-56*x^3*y^3*z^6*R+36*x^2*y*z^2+3*y*z^2-9*x^2*z^2-9*y^2*z^2+3*x^2*y-y*z+3*x^2*z-9*x^2*y^2+3*y^2*z)/R
	0
	,
--	3*y*y+ -- dy p term
--   (18*x*y*z^2-24*x^3*y*z^3+36*x^3*y*z^2+36*x^2*y*z^3-36*x^2*y^2*z+18*x*y^2*z+24*x^3*y^2*z^3-36*x^3*y^2*z^2-36*x^2*y^2*z^3+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+24*x^3*y^2*z+36*x^2*y^3*z-6*x^2*z^3+6*y^3*z-24*x*y^3*z-6*y^2*z^3-18*y^3*z^2+12*y^3*z^3+24*x*y^2*z^3+36*x*y^3*z^2-24*x*y^3*z^3-12*x*y*z^3+2*x*z^3-1620*x^4*y^5*z^5*R+1638*x^4*y^5*z^4*R+12*x*y^4*z+12*x^3*y^4*z-18*x^2*y^4*z+12*x*y^4*z^3-18*x*y^4*z^2-3*y^4*z+120*x^3*y^3*z^5*R+9*y^4*z^2-6*y^4*z^3-180*x^4*y^3*z^5*R+900*x^4*y^4*z^5*R+54*x^2*y^2*z^2-4*x^3*y^3*z^2*R-420*x^4*y^6*z^6*R+182*x^4*y^3*z^4*R-48*x^6*y^3*z^5*R+48*x^6*y^3*z^4*R+144*x^5*y^3*z^5*R-240*x^5*y^4*z^3*R+320*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+432*x^6*y^5*z^4*R+240*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-1296*x^5*y^5*z^4*R-720*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-144*x^6*y^5*z^3*R+80*x^6*y^4*z^3*R+432*x^5*y^5*z^3*R+432*x^3*y^5*z^3*R-576*x^4*y^5*z^3*R-600*x^3*y^4*z^5*R-1116*x^3*y^5*z^4*R+2*x^3*z+4*x^3*z^3-6*x^3*z^2-18*x^2*y^3+12*x^3*y^3-6*x^3*y^2+1080*x^3*y^5*z^5*R+6*y^3*x-3*y^4*x+9*y^4*x^2-6*y^4*x^3-1008*x^5*y^6*z^5*R+336*x^6*y^6*z^5*R-336*x^6*y^6*z^4*R+448*x^4*y^6*z^3*R+112*x^6*y^6*z^3*R-336*x^5*y^6*z^3*R+1260*x^4*y^6*z^5*R-1274*x^4*y^6*z^4*R+20*x^3*y^4*z^2*R-336*x^3*y^6*z^3*R+868*x^3*y^6*z^4*R-840*x^3*y^6*z^5*R+1008*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-36*x^3*y^5*z^2*R+240*x^5*y^4*z^6*R-432*x^5*y^5*z^6*R-80*x^6*y^4*z^6*R+144*x^6*y^5*z^6*R-48*x^5*y^3*z^6*R+16*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+28*x^3*y^6*z^2*R+280*x^3*y^6*z^6*R+18*x^4*y^5*z^2*R-14*x^4*y^6*z^2*R-144*x^5*y^3*z^4*R-124*x^3*y^3*z^4*R+620*x^3*y^4*z^4*R+2*x^4*y^3*z^2*R+48*x^3*y^3*z^3*R-10*x^4*y^4*z^2*R-240*x^3*y^4*z^3*R-300*x^4*y^4*z^6*R+540*x^4*y^5*z^6*R+60*x^4*y^3*z^6*R+200*x^3*y^4*z^6*R-360*x^3*y^5*z^6*R-40*x^3*y^3*z^6*R-54*x^2*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*z-3*x^2*z-3*x*y^2+9*x^2*y^2-3*y^2*z-10*x^2*y^4*z^2*R-16*x^2*y^3*z^3*R-80*x^3*y^7*z^6*R+2*x^2*y^3*z^2*R+112*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-266*x^2*y^6*z^4*R+252*x^2*y^6*z^5*R+180*x^2*y^4*z^5*R+364*x^4*y^7*z^4*R-360*x^4*y^7*z^5*R+288*x^5*y^7*z^5*R-288*x^5*y^7*z^4*R+96*x^3*y^7*z^3*R+96*x^5*y^7*z^3*R-128*x^4*y^7*z^3*R+240*x^3*y^7*z^5*R-248*x^3*y^7*z^4*R+18*x^2*y^5*z^2*R+38*x^2*y^3*z^4*R-32*x^2*y^7*z^3*R+76*x^2*y^7*z^4*R-72*x^2*y^7*z^5*R-14*x^2*y^6*z^2*R+120*x^4*y^7*z^6*R-96*x^5*y^7*z^6*R+4*x^2*y^7*z^2*R+24*x^2*y^7*z^6*R-8*x^3*y^7*z^2*R+4*x^4*y^7*z^2*R-96*x^6*y^7*z^5*R+96*x^6*y^7*z^4*R-32*x^6*y^7*z^3*R+32*x^6*y^7*z^6*R-190*x^2*y^4*z^4*R+342*x^2*y^5*z^4*R-36*x^2*y^3*z^5*R+80*x^2*y^4*z^3*R-144*x^2*y^5*z^3*R+12*x^2*y^3*z^6*R+108*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-60*x^2*y^4*z^6*R)/R
	0
   	,
--  3*z*z + -- dz p term
--   (18*x*y*z^2-24*x^3*y*z^3+24*x^3*y*z^2+36*x^2*y*z^3-54*x^2*y^2*z+18*x*y^2*z+24*x^3*y^3*z^2-36*x^3*y^2*z^2-36*x^2*y^3*z^2+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+36*x^3*y^2*z+36*x^2*y^3*z-18*x^2*z^3-12*x*y^3*z-18*y^2*z^3-6*y^3*z^2+12*y^3*z^3+36*x*y^2*z^3+24*x*y^3*z^2-24*x*y^3*z^3+6*y*z^3-24*x*y*z^3+6*x*z^3-1620*x^4*y^5*z^5*R+900*x^4*y^5*z^4*R+432*x^3*y^3*z^5*R-576*x^4*y^3*z^5*R+1638*x^4*y^4*z^5*R+54*x^2*y^2*z^2+9*x^2*z^4-3*x*z^4-4*x^3*y^2*z^3*R-420*x^4*y^6*z^6*R+320*x^4*y^3*z^4*R-144*x^6*y^3*z^5*R+80*x^6*y^3*z^4*R+432*x^5*y^3*z^5*R-144*x^5*y^4*z^3*R+182*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+240*x^6*y^5*z^4*R+432*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-720*x^5*y^5*z^4*R-1296*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-48*x^6*y^5*z^3*R+48*x^6*y^4*z^3*R+144*x^5*y^5*z^3*R+120*x^3*y^5*z^3*R-180*x^4*y^5*z^3*R-1116*x^3*y^4*z^5*R-600*x^3*y^5*z^4*R+12*x^3*z^3-6*x^3*z^2-6*x^2*y^3+4*x^3*y^3-6*x^3*y^2+2*x^3*y+1080*x^3*y^5*z^5*R-18*x^2*y*z^4+12*x^3*y*z^4+2*y^3*x-6*z^4*x^3+18*x^4*y^2*z^5*R-432*x^5*y^6*z^5*R+144*x^6*y^6*z^5*R-80*x^6*y^6*z^4*R+60*x^4*y^6*z^3*R+16*x^6*y^6*z^3*R-48*x^5*y^6*z^3*R+540*x^4*y^6*z^5*R-300*x^4*y^6*z^4*R+20*x^3*y^2*z^4*R-40*x^3*y^6*z^3*R+200*x^3*y^6*z^4*R-360*x^3*y^6*z^5*R+240*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-14*x^4*y^2*z^6*R+1008*x^5*y^4*z^6*R-1008*x^5*y^5*z^6*R-336*x^6*y^4*z^6*R+336*x^6*y^5*z^6*R-336*x^5*y^3*z^6*R+112*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+280*x^3*y^6*z^6*R-240*x^5*y^3*z^4*R-240*x^3*y^3*z^4*R-10*x^4*y^2*z^4*R+620*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+48*x^3*y^3*z^3*R+2*x^4*y^2*z^3*R-124*x^3*y^4*z^3*R+28*x^3*y^2*z^6*R-1274*x^4*y^4*z^6*R+1260*x^4*y^5*z^6*R+448*x^4*y^3*z^6*R+868*x^3*y^4*z^6*R-840*x^3*y^5*z^6*R-336*x^3*y^3*z^6*R-36*x^2*y*z^2-3*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*y-3*x^2*y-3*x*y^2+9*x^2*y^2-18*x*y^2*z^4+12*x*y^3*z^4+12*x*y*z^4+9*y^2*z^4-6*y^3*z^4-3*y*z^4-16*x^2*y^3*z^3*R+12*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-60*x^2*y^6*z^4*R+108*x^2*y^6*z^5*R+342*x^2*y^4*z^5*R+80*x^2*y^3*z^4*R-190*x^2*y^4*z^4*R+180*x^2*y^5*z^4*R-144*x^2*y^3*z^5*R+38*x^2*y^4*z^3*R-36*x^2*y^5*z^3*R+112*x^2*y^3*z^6*R+252*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-266*x^2*y^4*z^6*R+2*x^2*y^2*z^3*R-10*x^2*y^2*z^4*R-80*x^3*y^6*z^7*R+18*x^2*y^2*z^5*R+4*x^4*y^2*z^7*R-8*x^3*y^2*z^7*R+364*x^4*y^4*z^7*R-360*x^4*y^5*z^7*R-288*x^5*y^4*z^7*R+288*x^5*y^5*z^7*R-128*x^4*y^3*z^7*R+96*x^5*y^3*z^7*R+120*x^4*y^6*z^7*R-96*x^5*y^6*z^7*R+24*x^2*y^6*z^7*R+96*x^6*y^4*z^7*R-96*x^6*y^5*z^7*R-32*x^6*y^3*z^7*R+32*x^6*y^6*z^7*R-14*x^2*y^2*z^6*R+4*x^2*y^2*z^7*R-248*x^3*y^4*z^7*R+240*x^3*y^5*z^7*R+96*x^3*y^3*z^7*R+76*x^2*y^4*z^7*R-72*x^2*y^5*z^7*R-32*x^2*y^3*z^7*R)/R	
	0
end

if dim==2 then
	rhs = LuaUserVector("source2d")
else
	rhs = LuaUserVector("source3d")
end

NavierStokesDisc:set_source(rhs)

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)

op = AssembledOperator(domainDisc)

u = GridFunction(approxSpace)
u:set(0)

function StartValue_u2d(x,y,t) return y end
function StartValue_v2d(x,y,t) return -x end
function StartValue_p2d(x,y,t) return 0 end

function StartValue_u3d(x,y,z,t) return y end
function StartValue_v3d(x,y,z,t) return -x end
function StartValue_w3d(x,y,z,t) return 0 end
function StartValue_p3d(x,y,z,t) return 0 end

if dim==2 then
	Interpolate("StartValue_u2d", u, "u")
	Interpolate("StartValue_v2d", u, "v")
	Interpolate("StartValue_p2d", u, "p")
else
	Interpolate("StartValue_u3d", u, "u")
	Interpolate("StartValue_v3d", u, "v")
	Interpolate("StartValue_w3d", u, "w")
	Interpolate("StartValue_p3d", u, "p")
end

vanka = Vanka()
vanka:set_damp(1)

-- vanka = DiagVanka()

vankaSolver = LinearSolver()
vankaSolver:set_preconditioner(vanka)
vankaSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-12, true))

egsSolver = LinearSolver()
egsSolver:set_preconditioner(ElementGaussSeidel())
egsSolver:set_convergence_check(ConvCheck(10000, 1e-9, 1e-12, true))

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-7)
baseConvCheck:set_reduction(1e-3)
baseConvCheck:set_verbose(false)

vankaBase = LinearSolver()
vankaBase:set_preconditioner(DiagVanka())
vankaBase:set_convergence_check(baseConvCheck)

crilut = CRILUT()
crilut:set_threshold(1e-1,1e-4,1e-4,1e-4)
crilut:set_damp(1)
-- crilut:set_info(true)
ilutBase = ILUT()
ilutBase:set_threshold(1e-4)
ilutBase:set_info(false)
ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(ilutBase)
ilutSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-12, false))


gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(vankaBase)
if discType=="fv1" then
	ilut = ILUT()
	ilut:set_threshold(1e-3)
	ilut:set_info(false)
	gmg:set_smoother(ilut)
end
if discType=="fvcr" then
	gmg:set_smoother(vanka)
end
gmg:set_cycle_type(1)
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
gmg:set_damp(MinimalResiduumDamping())
-- gmg:add_prolongation_post_process(AverageComponent("p"))
-- gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)
-- create Linear Solver
BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
BiCGStabSolver:set_convergence_check(ConvCheck(100000, 1e-7, 1e-12, true))

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(ConvCheck(10000, 1e-10, 1e-12, true))


-- choose a solver
if discType=="fvcr" or discType=="fv1" then
	solver = gmgSolver
end
if discType=="fe" or discType=="fv" then
	solver = egsSolver
end

newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(100)
newtonConvCheck:set_minimum_defect(1e-7)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(30)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.85)
newtonLineSearch:set_accept_best(true)

dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

tBefore = os.clock()

newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end

SaveVectorForConnectionViewer(u, "StartSolution.vec")

if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed.") exit()
end

l2error = L2Error(uSolution, u, "u", 0.0, 1, "Inner")
write("L2Error in u component is "..l2error .."\n")
l2error = L2Error(vSolution, u, "v", 0.0, 1, "Inner")
write("L2Error in v component is "..l2error .."\n")
l2error = L2Error(pSolution, u, "p", 0.0, 1, "Inner")
write("L2Error in p component is "..l2error .."\n")


out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_element("u,v", "velocity")
out:select_element("u", "u")
out:select_element("v", "v")
out:select_element("p", "p")
out:print("StokesSolution", u)

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.")
