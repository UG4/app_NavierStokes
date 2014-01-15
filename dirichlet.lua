--------------------------------------------------------------------------------
--
--  Lua - Script to test the Navier-Stokes implementation
--
--  Author: Christian Wehner
--
--	A theoretical example to test the Navier-Stokes discretization.
--	The boundary conditions are inflow boundary conditions 
--  (Dirichlet conditions for the velocity) on the whole boundary.
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
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim 		= util.GetParamNumber("-dim", 2)
type 		= util.GetParam("-type", "fvcr")
order 		= util.GetParamNumber("-order", 1)
vorder 		= util.GetParamNumber("-vorder", order)
porder 		= util.GetParamNumber("-porder", vorder-1)
R 			= util.GetParamNumber("-R", 1)
bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "no", "Upwind type")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "COR", "Diffusion length type")
linred      = util.GetParam("-linred", 1e-2 , "Linear reduction")
nlintol     = util.GetParam("-nlintol", 1e-8, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)
numRefs 	= util.GetParamNumber("-numRefs",2)

if 	dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_1x1.ugx")
else
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_tets.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_1x1x1.ugx")
end
if dim~=2 and dim~=3 then
   print("Chosen Dimension " .. dim .. " not supported. Exiting.") exit() 
end

InitUG(dim, AlgebraType("CPU", 1))



print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    type       = " .. type)
print("    grid       	= " .. gridName)
print("    v ansatz order = " ..vorder)
print("    p ansatz order = " ..porder)
print("    no laplace       = " .. tostring(bNoLaplace))
print("    exact jacobian   = " .. tostring(bExactJac))
print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    pac upwind       = " .. tostring(bPac))
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)
print("    linear reduction = " .. linred)
print("    linear tolerance = " .. lintol)
print("    nonlinear reduction = " .. nlinred)
print("    nonlinear tolerance = " .. nlintol)
print("    Reynolds number = " .. R)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

requiredSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
approxSpace = ApproximationSpace(dom)

-- components of the velocity 
if 		dim == 1 then VelCmp = {"u"} FctCmp = {"u", "p"}
elseif  dim == 2 then VelCmp = {"u", "v"} FctCmp = {"u", "v", "p"}
elseif  dim == 3 then VelCmp = {"u", "v", "w"} FctCmp = {"u", "v", "w", "p"}
else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end

if type == "fv1" then
	approxSpace:add_fct(FctCmp, "Lagrange", 1) 
elseif type == "fv" then
	approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
	if porder==0 then
		print("p order 0 not possible for fv type. Exiting.") exit()
	else
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
elseif type == "fe" then
	if porder==0 then
		if vorder==1 then
			approxSpace:add_fct(VelCmp, "Crouzeix-Raviart",1)
		else
			approxSpace:add_fct(VelCmp, "Lagrange", vorder)
		end
		approxSpace:add_fct("p", "piecewise-constant") 
	else
		approxSpace:add_fct(VelCmp, "Lagrange", vorder) 
		approxSpace:add_fct("p", "Lagrange", porder) 
	end
elseif type=="fvcr" then
	approxSpace:add_fct(VelCmp, "Crouzeix-Raviart")
	approxSpace:add_fct("p", "piecewise-constant") 
else print("Disc Type '"..type.."' not supported."); exit(); end

-- print statistic on the distributed dofs
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

NavierStokesDisc = NavierStokes(fctUsed, "Inner", type)
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_laplace( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity(1.0/R);

--upwind if available
if type == "fv1" or type == "fvcr" then
	NavierStokesDisc:set_upwind(upwind)
	NavierStokesDisc:set_peclet_blend(bPecletBlend)
end

-- fv1 must be stablilized
if type == "fv1" then
	NavierStokesDisc:set_stabilization(stab, diffLength)
	NavierStokesDisc:set_pac_upwind(bPac)
end

-- fe must be stabilized for (Pk, Pk) space
if type == "fe" and porder == vorder then
	NavierStokesDisc:set_stabilization(3)
end

--------------------------------------------------------------------------------
-- Boundary conditions
--------------------------------------------------------------------------------

function usol2d(x, y, t)
	return 2*x^2*(1-x)^2*(y*(1-y)^2-y^2*(1-y))
--	return y
end

function vsol2d(x,y,t)
	return -2*y^2*(1-y)^2*(x*(1-x)^2-x^2*(1-x))
--	return -x
end

function psol2d(x,y,t)
	return x^2+y^2
--	return 0
--	return x+y
end

function usol3d(x,y,z,t)
	return 2*x^2*y*z*(2*z-1)*(z-1)*(2*y-1)*(y-1)*(x-1)^2
--	return y
end

function vsol3d(x,y,z,t)
	return -x*y^2*z*(2*z-1)*(z-1)*(y-1)^2*(2*x-1)*(x-1)
--	return -x
end

function wsol3d(x,y,z,t)
	return -x*y*z^2*(z-1)^2*(2*y-1)*(y-1)*(2*x-1)*(x-1)
--	return 0
end

function psol3d(x,y,z,t)
	return x^2+y^2+z^2
--	return 0
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t),vsol2d(x, y, t)
end

function inletVel3d(x, y,z, t)
	return usol3d(x,y,z,t),vsol3d(x,y,z,t), wsol3d(x,y,z,t)
end

InletDisc = NavierStokesInflow(NavierStokesDisc)
InletDisc:add("inletVel"..dim.."d", "Boundary")

--------------------------------------------------------------------------------
-- Source
--------------------------------------------------------------------------------

if bStokes == true then
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
else
	function source2d(x, y, t)
		return 
		2*x+ -- dx p term
	    (-4*y+12*x^2+12*y^2-24*x^3-8*y^3+24*x*y-72*x*y^2+48*x*y^3-48*x^2*y+72*x^2*y^2-48*x^2*y^3+48*x^3*y-24*x^4*y+12*x^4-16*x^3*y^3*R-20*x^4*y^2*R+36*x^5*y^2*R+28*x^3*y^4*R-28*x^6*y^2*R-24*x^3*y^5*R+80*x^4*y^3*R-140*x^4*y^4*R+120*x^4*y^5*R-144*x^5*y^3*R+252*x^5*y^4*R-216*x^5*y^5*R+112*x^6*y^3*R-32*x^7*y^3*R+8*x^7*y^2*R-196*x^6*y^4*R+168*x^6*y^5*R-40*x^4*y^6*R+72*x^5*y^6*R-56*x^6*y^6*R+56*x^7*y^4*R-48*x^7*y^5*R+16*x^7*y^6*R+8*x^3*y^6*R+4*x^3*y^2*R)/R
	    --	-x
		,
		2*y+ -- dy p term
		(4*x-12*x^2-12*y^2+8*x^3+24*y^3-24*x*y+48*x*y^2-48*x*y^3+72*x^2*y-72*x^2*y^2-48*x^3*y+48*x^3*y^2+24*x*y^4-12*y^4-16*x^3*y^3*R+80*x^3*y^4*R-144*x^3*y^5*R+28*x^4*y^3*R-140*x^4*y^4*R+252*x^4*y^5*R-24*x^5*y^3*R+120*x^5*y^4*R-216*x^5*y^5*R+8*x^6*y^3*R-40*x^6*y^4*R+72*x^6*y^5*R-196*x^4*y^6*R+168*x^5*y^6*R-56*x^6*y^6*R+112*x^3*y^6*R-20*x^2*y^4*R+36*x^2*y^5*R-28*x^2*y^6*R-32*x^3*y^7*R+56*x^4*y^7*R-48*x^5*y^7*R+16*x^6*y^7*R+8*x^2*y^7*R+4*x^2*y^3*R)/R
		--	-y
	end	
	
	function source3d(x, y,z, t)
		return 
		2*x +  -- dx p term
	    2*(-18*x*y*z^2-12*x^4*y*z^3+18*x^4*y*z^2+24*x^3*y*z^3-36*x^3*y*z^2-24*x^2*y*z^3+36*x^2*y^2*z-18*x*y^2*z-24*x^2*y^3*z^3+36*x^2*y^3*z^2+36*x^2*y^2*z^3-18*x^2*y*z+6*x*y*z+24*x^3*y*z-12*x^4*y*z+54*x*y^2*z^2-12*x^4*y^3*z+18*x^4*y^2*z+24*x^3*y^3*z-36*x^3*y^2*z-24*x^2*y^3*z+6*x^2*z^3-2*y^3*z+12*x*y^3*z+6*y^2*z^3+6*y^3*z^2-4*y^3*z^3-36*x*y^2*z^3-36*x*y^3*z^2+24*x*y^3*z^3-2*y*z^3+12*x*y*z^3-1440*x^4*y^5*z^5*R+1620*x^4*y^5*z^4*R+168*x^3*y^3*z^5*R-840*x^4*y^3*z^5*R+1620*x^4*y^4*z^5*R-54*x^2*y^2*z^2+36*x^5*y^2*z^2*R-20*x^4*y^2*z^2*R-20*x^3*y^3*z^2*R-20*x^3*y^2*z^3*R-160*x^4*y^6*z^6*R+940*x^4*y^3*z^4*R+336*x^7*y^3*z^5*R-376*x^7*y^3*z^4*R-1176*x^6*y^3*z^5*R+1316*x^6*y^3*z^4*R+1512*x^5*y^3*z^5*R-1692*x^5*y^4*z^3*R+940*x^4*y^4*z^3*R+576*x^7*y^5*z^5*R-648*x^7*y^5*z^4*R-648*x^7*y^4*z^5*R+728*x^7*y^4*z^4*R-2016*x^6*y^5*z^5*R+2268*x^6*y^5*z^4*R+2268*x^6*y^4*z^5*R-2548*x^6*y^4*z^4*R+2592*x^5*y^5*z^5*R-2916*x^5*y^5*z^4*R-2916*x^5*y^4*z^5*R+864*x^5*y^3*z^3*R-480*x^4*y^3*z^3*R-672*x^6*y^3*z^3*R+192*x^7*y^3*z^3*R-1820*x^4*y^4*z^4*R+336*x^7*y^5*z^3*R-376*x^7*y^4*z^3*R-1176*x^6*y^5*z^3*R+1316*x^6*y^4*z^3*R+1512*x^5*y^5*z^3*R-324*x^5*y^2*z^5*R+168*x^3*y^5*z^3*R-840*x^4*y^5*z^3*R-324*x^3*y^4*z^5*R-324*x^3*y^5*z^4*R-6*x^3*z+3*x^4*z-12*x^3*z^3+6*x^4*z^3+18*x^3*z^2-9*x^4*z^2+6*x^2*y^3-12*x^3*y^3+6*x^4*y^3+18*x^3*y^2-9*x^4*y^2-6*x^3*y+3*x^4*y+288*x^3*y^5*z^5*R+4*x^3*y^2*z^2*R+180*x^4*y^2*z^5*R-864*x^5*y^6*z^5*R+672*x^6*y^6*z^5*R-756*x^6*y^6*z^4*R+280*x^4*y^6*z^3*R+392*x^6*y^6*z^3*R-504*x^5*y^6*z^3*R+480*x^4*y^6*z^5*R-540*x^4*y^6*z^4*R+40*x^3*y^4*z^2*R+40*x^3*y^2*z^4*R-56*x^3*y^6*z^3*R+108*x^3*y^6*z^4*R-96*x^3*y^6*z^5*R+972*x^5*y^6*z^4*R+3276*x^5*y^4*z^4*R+108*x^5*y^2*z^6*R-60*x^4*y^2*z^6*R+140*x^6*y^2*z^3*R-40*x^7*y^2*z^3*R+252*x^6*y^2*z^5*R-72*x^7*y^2*z^5*R-280*x^6*y^2*z^4*R+80*x^7*y^2*z^4*R-324*x^5*y^5*z^2*R+252*x^6*y^5*z^2*R-72*x^7*y^5*z^2*R-280*x^6*y^4*z^2*R+80*x^7*y^4*z^2*R+140*x^6*y^3*z^2*R-40*x^7*y^3*z^2*R-36*x^3*y^5*z^2*R+972*x^5*y^4*z^6*R-864*x^5*y^5*z^6*R-756*x^6*y^4*z^6*R+672*x^6*y^5*z^6*R-504*x^5*y^3*z^6*R+392*x^6*y^3*z^6*R+288*x^5*y^6*z^6*R-224*x^6*y^6*z^6*R+12*x^3*y^6*z^2*R-28*x^6*y^2*z^2*R+32*x^3*y^6*z^6*R+180*x^4*y^5*z^2*R-60*x^4*y^6*z^2*R+108*x^5*y^6*z^2*R-84*x^6*y^6*z^2*R-84*x^6*y^2*z^6*R+8*x^7*y^2*z^2*R-192*x^7*y^6*z^5*R+216*x^7*y^6*z^4*R-112*x^7*y^6*z^3*R+216*x^7*y^4*z^6*R-192*x^7*y^5*z^6*R-112*x^7*y^3*z^6*R+64*x^7*y^6*z^6*R+24*x^7*y^6*z^2*R+24*x^7*y^2*z^6*R-1692*x^5*y^3*z^4*R-188*x^3*y^3*z^4*R-200*x^4*y^2*z^4*R+360*x^5*y^2*z^4*R+364*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+100*x^4*y^3*z^2*R-180*x^5*y^3*z^2*R+96*x^3*y^3*z^3*R+100*x^4*y^2*z^3*R-180*x^5*y^2*z^3*R-200*x^4*y^4*z^2*R+360*x^5*y^4*z^2*R-188*x^3*y^4*z^3*R+12*x^3*y^2*z^6*R-540*x^4*y^4*z^6*R+480*x^4*y^5*z^6*R+280*x^4*y^3*z^6*R+108*x^3*y^4*z^6*R-96*x^3*y^5*z^6*R-56*x^3*y^3*z^6*R+36*x^2*y*z^2+3*y*z^2-9*x^2*z^2-9*y^2*z^2+3*x^2*y-y*z+3*x^2*z-9*x^2*y^2+3*y^2*z)/R
	--	-x
		,
		2*y + -- dy p term
	   (18*x*y*z^2-24*x^3*y*z^3+36*x^3*y*z^2+36*x^2*y*z^3-36*x^2*y^2*z+18*x*y^2*z+24*x^3*y^2*z^3-36*x^3*y^2*z^2-36*x^2*y^2*z^3+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+24*x^3*y^2*z+36*x^2*y^3*z-6*x^2*z^3+6*y^3*z-24*x*y^3*z-6*y^2*z^3-18*y^3*z^2+12*y^3*z^3+24*x*y^2*z^3+36*x*y^3*z^2-24*x*y^3*z^3-12*x*y*z^3+2*x*z^3-1620*x^4*y^5*z^5*R+1638*x^4*y^5*z^4*R+12*x*y^4*z+12*x^3*y^4*z-18*x^2*y^4*z+12*x*y^4*z^3-18*x*y^4*z^2-3*y^4*z+120*x^3*y^3*z^5*R+9*y^4*z^2-6*y^4*z^3-180*x^4*y^3*z^5*R+900*x^4*y^4*z^5*R+54*x^2*y^2*z^2-4*x^3*y^3*z^2*R-420*x^4*y^6*z^6*R+182*x^4*y^3*z^4*R-48*x^6*y^3*z^5*R+48*x^6*y^3*z^4*R+144*x^5*y^3*z^5*R-240*x^5*y^4*z^3*R+320*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+432*x^6*y^5*z^4*R+240*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-1296*x^5*y^5*z^4*R-720*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-144*x^6*y^5*z^3*R+80*x^6*y^4*z^3*R+432*x^5*y^5*z^3*R+432*x^3*y^5*z^3*R-576*x^4*y^5*z^3*R-600*x^3*y^4*z^5*R-1116*x^3*y^5*z^4*R+2*x^3*z+4*x^3*z^3-6*x^3*z^2-18*x^2*y^3+12*x^3*y^3-6*x^3*y^2+1080*x^3*y^5*z^5*R+6*y^3*x-3*y^4*x+9*y^4*x^2-6*y^4*x^3-1008*x^5*y^6*z^5*R+336*x^6*y^6*z^5*R-336*x^6*y^6*z^4*R+448*x^4*y^6*z^3*R+112*x^6*y^6*z^3*R-336*x^5*y^6*z^3*R+1260*x^4*y^6*z^5*R-1274*x^4*y^6*z^4*R+20*x^3*y^4*z^2*R-336*x^3*y^6*z^3*R+868*x^3*y^6*z^4*R-840*x^3*y^6*z^5*R+1008*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-36*x^3*y^5*z^2*R+240*x^5*y^4*z^6*R-432*x^5*y^5*z^6*R-80*x^6*y^4*z^6*R+144*x^6*y^5*z^6*R-48*x^5*y^3*z^6*R+16*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+28*x^3*y^6*z^2*R+280*x^3*y^6*z^6*R+18*x^4*y^5*z^2*R-14*x^4*y^6*z^2*R-144*x^5*y^3*z^4*R-124*x^3*y^3*z^4*R+620*x^3*y^4*z^4*R+2*x^4*y^3*z^2*R+48*x^3*y^3*z^3*R-10*x^4*y^4*z^2*R-240*x^3*y^4*z^3*R-300*x^4*y^4*z^6*R+540*x^4*y^5*z^6*R+60*x^4*y^3*z^6*R+200*x^3*y^4*z^6*R-360*x^3*y^5*z^6*R-40*x^3*y^3*z^6*R-54*x^2*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*z-3*x^2*z-3*x*y^2+9*x^2*y^2-3*y^2*z-10*x^2*y^4*z^2*R-16*x^2*y^3*z^3*R-80*x^3*y^7*z^6*R+2*x^2*y^3*z^2*R+112*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-266*x^2*y^6*z^4*R+252*x^2*y^6*z^5*R+180*x^2*y^4*z^5*R+364*x^4*y^7*z^4*R-360*x^4*y^7*z^5*R+288*x^5*y^7*z^5*R-288*x^5*y^7*z^4*R+96*x^3*y^7*z^3*R+96*x^5*y^7*z^3*R-128*x^4*y^7*z^3*R+240*x^3*y^7*z^5*R-248*x^3*y^7*z^4*R+18*x^2*y^5*z^2*R+38*x^2*y^3*z^4*R-32*x^2*y^7*z^3*R+76*x^2*y^7*z^4*R-72*x^2*y^7*z^5*R-14*x^2*y^6*z^2*R+120*x^4*y^7*z^6*R-96*x^5*y^7*z^6*R+4*x^2*y^7*z^2*R+24*x^2*y^7*z^6*R-8*x^3*y^7*z^2*R+4*x^4*y^7*z^2*R-96*x^6*y^7*z^5*R+96*x^6*y^7*z^4*R-32*x^6*y^7*z^3*R+32*x^6*y^7*z^6*R-190*x^2*y^4*z^4*R+342*x^2*y^5*z^4*R-36*x^2*y^3*z^5*R+80*x^2*y^4*z^3*R-144*x^2*y^5*z^3*R+12*x^2*y^3*z^6*R+108*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-60*x^2*y^4*z^6*R)/R
	--	-y
	   	,
	     2*z + -- dz p term
	    (18*x*y*z^2-24*x^3*y*z^3+24*x^3*y*z^2+36*x^2*y*z^3-54*x^2*y^2*z+18*x*y^2*z+24*x^3*y^3*z^2-36*x^3*y^2*z^2-36*x^2*y^3*z^2+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+36*x^3*y^2*z+36*x^2*y^3*z-18*x^2*z^3-12*x*y^3*z-18*y^2*z^3-6*y^3*z^2+12*y^3*z^3+36*x*y^2*z^3+24*x*y^3*z^2-24*x*y^3*z^3+6*y*z^3-24*x*y*z^3+6*x*z^3-1620*x^4*y^5*z^5*R+900*x^4*y^5*z^4*R+432*x^3*y^3*z^5*R-576*x^4*y^3*z^5*R+1638*x^4*y^4*z^5*R+54*x^2*y^2*z^2+9*x^2*z^4-3*x*z^4-4*x^3*y^2*z^3*R-420*x^4*y^6*z^6*R+320*x^4*y^3*z^4*R-144*x^6*y^3*z^5*R+80*x^6*y^3*z^4*R+432*x^5*y^3*z^5*R-144*x^5*y^4*z^3*R+182*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+240*x^6*y^5*z^4*R+432*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-720*x^5*y^5*z^4*R-1296*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-48*x^6*y^5*z^3*R+48*x^6*y^4*z^3*R+144*x^5*y^5*z^3*R+120*x^3*y^5*z^3*R-180*x^4*y^5*z^3*R-1116*x^3*y^4*z^5*R-600*x^3*y^5*z^4*R+12*x^3*z^3-6*x^3*z^2-6*x^2*y^3+4*x^3*y^3-6*x^3*y^2+2*x^3*y+1080*x^3*y^5*z^5*R-18*x^2*y*z^4+12*x^3*y*z^4+2*y^3*x-6*z^4*x^3+18*x^4*y^2*z^5*R-432*x^5*y^6*z^5*R+144*x^6*y^6*z^5*R-80*x^6*y^6*z^4*R+60*x^4*y^6*z^3*R+16*x^6*y^6*z^3*R-48*x^5*y^6*z^3*R+540*x^4*y^6*z^5*R-300*x^4*y^6*z^4*R+20*x^3*y^2*z^4*R-40*x^3*y^6*z^3*R+200*x^3*y^6*z^4*R-360*x^3*y^6*z^5*R+240*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-14*x^4*y^2*z^6*R+1008*x^5*y^4*z^6*R-1008*x^5*y^5*z^6*R-336*x^6*y^4*z^6*R+336*x^6*y^5*z^6*R-336*x^5*y^3*z^6*R+112*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+280*x^3*y^6*z^6*R-240*x^5*y^3*z^4*R-240*x^3*y^3*z^4*R-10*x^4*y^2*z^4*R+620*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+48*x^3*y^3*z^3*R+2*x^4*y^2*z^3*R-124*x^3*y^4*z^3*R+28*x^3*y^2*z^6*R-1274*x^4*y^4*z^6*R+1260*x^4*y^5*z^6*R+448*x^4*y^3*z^6*R+868*x^3*y^4*z^6*R-840*x^3*y^5*z^6*R-336*x^3*y^3*z^6*R-36*x^2*y*z^2-3*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*y-3*x^2*y-3*x*y^2+9*x^2*y^2-18*x*y^2*z^4+12*x*y^3*z^4+12*x*y*z^4+9*y^2*z^4-6*y^3*z^4-3*y*z^4-16*x^2*y^3*z^3*R+12*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-60*x^2*y^6*z^4*R+108*x^2*y^6*z^5*R+342*x^2*y^4*z^5*R+80*x^2*y^3*z^4*R-190*x^2*y^4*z^4*R+180*x^2*y^5*z^4*R-144*x^2*y^3*z^5*R+38*x^2*y^4*z^3*R-36*x^2*y^5*z^3*R+112*x^2*y^3*z^6*R+252*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-266*x^2*y^4*z^6*R+2*x^2*y^2*z^3*R-10*x^2*y^2*z^4*R-80*x^3*y^6*z^7*R+18*x^2*y^2*z^5*R+4*x^4*y^2*z^7*R-8*x^3*y^2*z^7*R+364*x^4*y^4*z^7*R-360*x^4*y^5*z^7*R-288*x^5*y^4*z^7*R+288*x^5*y^5*z^7*R-128*x^4*y^3*z^7*R+96*x^5*y^3*z^7*R+120*x^4*y^6*z^7*R-96*x^5*y^6*z^7*R+24*x^2*y^6*z^7*R+96*x^6*y^4*z^7*R-96*x^6*y^5*z^7*R-32*x^6*y^3*z^7*R+32*x^6*y^6*z^7*R-14*x^2*y^2*z^6*R+4*x^2*y^2*z^7*R-248*x^3*y^4*z^7*R+240*x^3*y^5*z^7*R+96*x^3*y^3*z^7*R+76*x^2*y^4*z^7*R-72*x^2*y^5*z^7*R-32*x^2*y^3*z^7*R)/R	
	--	0
	end	
end

NavierStokesDisc:set_source("source"..dim.."d")

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(InletDisc)

op = AssembledOperator(domainDisc)

u = GridFunction(approxSpace)
u:set(0)

vanka = Vanka()
vanka:set_damp(1)

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

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(ConvCheck(100, nlintol, nlinred, true))
newtonSolver:set_line_search(StandardLineSearch(30, 1.0, 0.85, true))
newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

tBefore = os.clock()

newtonSolver:init(op)

if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end

SaveVectorForConnectionViewer(u, "StartSolution.vec")

if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed.") exit()
end

l2error = L2Error("usol"..dim.."d", u, "u", 0.0, 1, "Inner")
write("L2Error in u component is "..l2error .."\n")
l2error = L2Error("vsol"..dim.."d", u, "v", 0.0, 1, "Inner")
write("L2Error in v component is "..l2error .."\n")
if dim==3 then
	l2error = L2Error("wsol"..dim.."d", u, "w", 0.0, 1, "Inner")
	write("L2Error in w component is "..l2error .."\n")
end
-- to make error computation for p reasonable
-- p would have to be adjusted by adding a reasonable constant
-- l2error = L2Error("psol"..dim.."d", u, "p", 0.0, 1, "Inner")
-- write("L2Error in p component is "..l2error .."\n")
maxerror = MaxError("usol"..dim.."d", u, "u")
write("Maximum error in u component is "..maxerror .."\n")
maxerror = MaxError("vsol"..dim.."d", u, "v")
write("Maximum error in v component is "..maxerror .."\n")
if dim==3 then
	maxerror = MaxError("wsol"..dim.."d", u, "w")
	write("Maximum error in w component is "..maxerror .."\n")
end

vtkWriter = VTKOutput()
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "pressure")
vtkWriter:print("Dirichlet", u)

tAfter = os.clock()
print("Computation took " .. tAfter-tBefore .. " seconds.")
