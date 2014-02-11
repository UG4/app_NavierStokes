--------------------------------------------------------------------------------
--
--  Lua - Script to test the Navier-Stokes implementation
--
--  Author: Christian Wehner
--
--	A theoretical example to test the Navier-Stokes discretization.
--	The boundary conditions are inflow boundary conditions 
--  (Dirichlet conditions for the velocity) on the whole boundary.
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

dim 		= util.GetParamNumber("-dim", 2)
numRefs 	= util.GetParamNumber("-numRefs",2)
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)
bConvRates  = util.HasParamOption("-convRate", "compute convergence rates")

bInstat     = util.HasParamOption("-instat", "time-dependent solution")
if bInstat then
dt 			= util.GetParamNumber("-dt", 0.05)
numTimeSteps= util.GetParamNumber("-numTimeSteps", 50)
end

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
nlintol     = util.GetParam("-nlintol", 1e-10, "Nonlinear tolerance")
lintol      = util.GetParam("-lintol", nlintol*0.5, "Linear tolerance")
nlinred     = util.GetParam("-nlinred", nlintol*0.1, "Nonlinear reduction")

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

discType, vorder, porder = util.ns.parseParams()

print(" Chosen Parameters:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)
print("    discType       = " .. discType)
if discType == "fv" or discType == "fe" then
print("    v ansatz order = " ..vorder)
print("    p ansatz order = " ..porder)
end
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
-- Source
--------------------------------------------------------------------------------
--[[
u:=sin(2*pi*(x+t))*cos(2*pi*y);
v:=-cos(2*pi*(x+t))*sin(2*pi*y);
# chose p
p:=0.25*x*x;
# rhs is chosen so that Navier-Stokes system is fulfilled
rhsu:=factor(diff(u,t)+simplify(-1/R*diff(diff(u,x),x)-1/R*diff(diff(u,y),y)+u*diff(u,x)+v*diff(u,y)+diff(p,x)));
rhsv:=factor(diff(v,t)+simplify(-1/R*diff(diff(v,x),x)-1/R*diff(diff(v,y),y)+u*diff(v,x)+v*diff(v,y)+diff(p,y)));
--]]
if bInstat then 

	function usol2d(x, y, t) return math.sin(2*math.pi*(x+t))*math.cos(2*math.pi*y)	end
	function vsol2d(x, y, t) return -math.cos(2*math.pi*(x+t))*math.sin(2*math.pi*y)end
	function psol2d(x, y, t) return math.sin(x)*math.sin(y)*math.sin(t) end

	function source2d(x, y, t)
		return 
		2*math.pi*(math.cos(2*math.pi*(x+t))*math.cos(2*math.pi*y)*R+4*math.sin(2*math.pi*(x+t))*math.pi*math.cos(2*math.pi*y)+math.sin(2*math.pi*(x+t))*math.cos(2*math.pi*(x+t))*R)/R
		,
		-2*math.pi*math.sin(2*math.pi*y)*(-math.sin(2*math.pi*(x+t))*R+4*math.cos(2*math.pi*(x+t))*math.pi-math.cos(2*math.pi*y)*R)/R
	end

--[[
The analytical solution can be constructed e.g. via maple:
restart:
with(codegen,C):
# choose divergence free velocity (using g) and pressure p 
g:=x^a+y^a;
u:=diff(g,y);
v:=-diff(g,x); 
p:=x^b+y^b-(2/(b+1));

# rhs is chosen so that (Navier)-Stokes system is fulfilled

laplace_u := factor(simplify(-1/R*(diff(u,x,x)+diff(u,y,y)))); 
nonlin_u := factor(simplify(u*diff(u,x)+v*diff(u,y))); 
press_u := factor(simplify(diff(p,x)));

laplace_v := factor(simplify(-1/R*(diff(v,x,x)+diff(v,y,y))));
nonlin_v := factor(simplify(u*diff(v,x)+v*diff(v,y)));
press_v := factor(simplify(diff(p,y)));

C( u ); 
C( v ); 
C( p ); 
C( laplace_u + nonlin_u + press_u );
C( laplace_v + nonlin_v + press_v );
 --]]
else
	PI = math.pi
	function usol2d(x, y, t) return 2*math.cos(2*PI*y)*PI  end
	function vsol2d(x, y, t) return -2*math.cos(2*PI*x)*PI end
	function psol2d(x, y, t) return math.sin(2*PI*x)+math.sin(2*PI*y)  end

	function ugrad2d(x, y, t) return 0, -4*PI*PI*math.sin(2*PI*y)  end
	function vgrad2d(x, y, t) return 4*PI*PI*math.sin(2*PI*x), 0 end
	function pgrad2d(x, y, t) return 2*PI*math.cos(2*PI*x),2*PI*math.cos(2*PI*y)  end

	if bStokes == true then
		function source2d(x, y, t)
			return 
			8/R*math.cos(2*PI*y)*PI*PI*PI+2*math.cos(2*PI*x)*PI,
			-8/R*math.cos(2*PI*x)*PI*PI*PI+2*math.cos(2*PI*y)*PI
		end
	else
		function source2d(x, y, t)
			return 
			8.0/R*math.cos(2.0*PI*y)*PI*PI*PI+8.0*math.cos(2.0*PI*x)*PI*PI*PI*math.sin(2.0*PI*y)+2.0*math.cos(2.0*PI*x)*PI,
 			-8.0/R*math.cos(2.0*PI*x)*PI*PI*PI+8.0*math.cos(2.0*PI*y)*PI*PI*PI*math.sin(2.0*PI*x)+2.0*math.cos(2.0*PI*y)*PI
		end	
	end

--[[
	a=vorder
	b=porder

	function usol2d(x, y, t) return math.pow(y,(a-1))*a  end
	function vsol2d(x, y, t) return -math.pow(x,(a-1))*a end
	function psol2d(x, y, t) return math.pow(x,b)+math.pow(y,b)-(2/(b+1))  end

	if bStokes == true then
		function source2d(x, y, t)
			return 
			-math.pow(y,a-3)*a*(a-1)*(a-2)/R+math.pow(x,b-1)*b,
			math.pow(x,a-3)*a*(a-1)*(a-2)/R+math.pow(y,b-1)*b
		end
	else
		function source2d(x, y, t)
			return 
			-math.pow(y,a-3)*a*(a-1)*(a-2)/R-a*a*math.pow(x,a-1)*math.pow(y,a-2)*(a-1)+math.pow(x,b-1)*b,
			math.pow(x,a-3)*a*(a-1)*(a-2)/R-a*a*math.pow(y,a-1)*math.pow(x,a-2)*(a-1)+math.pow(y,b-1)*b
		end	
	end
--]]

	function usol3d(x,y,z,t) return 2*x^2*y*z*(2*z-1)*(z-1)*(2*y-1)*(y-1)*(x-1)^2 end
	function vsol3d(x,y,z,t) return -x*y^2*z*(2*z-1)*(z-1)*(y-1)^2*(2*x-1)*(x-1)end
	function wsol3d(x,y,z,t) return -x*y*z^2*(z-1)^2*(2*y-1)*(y-1)*(2*x-1)*(x-1)end
	function psol3d(x,y,z,t) return x^2+y^2+z^2 end
	
	if bStokes == true then
		function source3d(x, y,z, t)
			return 
			3*x*x+ -- dx p term
			2*(-18*x*y*z^2-12*x^4*y*z^3+18*x^4*y*z^2+24*x^3*y*z^3-36*x^3*y*z^2-24*x^2*y*z^3+36*x^2*y^2*z-18*x*y^2*z-24*x^2*y^3*z^3+36*x^2*y^3*z^2+36*x^2*y^2*z^3-18*x^2*y*z+6*x*y*z+24*x^3*y*z-12*x^4*y*z+54*x*y^2*z^2-12*x^4*y^3*z+18*x^4*y^2*z+24*x^3*y^3*z-36*x^3*y^2*z-24*x^2*y^3*z+6*x^2*z^3-2*y^3*z+12*x*y^3*z+6*y^2*z^3+6*y^3*z^2-4*y^3*z^3-36*x*y^2*z^3-36*x*y^3*z^2+24*x*y^3*z^3-2*y*z^3+12*x*y*z^3-1440*x^4*y^5*z^5*R+1620*x^4*y^5*z^4*R+168*x^3*y^3*z^5*R-840*x^4*y^3*z^5*R+1620*x^4*y^4*z^5*R-54*x^2*y^2*z^2+36*x^5*y^2*z^2*R-20*x^4*y^2*z^2*R-20*x^3*y^3*z^2*R-20*x^3*y^2*z^3*R-160*x^4*y^6*z^6*R+940*x^4*y^3*z^4*R+336*x^7*y^3*z^5*R-376*x^7*y^3*z^4*R-1176*x^6*y^3*z^5*R+1316*x^6*y^3*z^4*R+1512*x^5*y^3*z^5*R-1692*x^5*y^4*z^3*R+940*x^4*y^4*z^3*R+576*x^7*y^5*z^5*R-648*x^7*y^5*z^4*R-648*x^7*y^4*z^5*R+728*x^7*y^4*z^4*R-2016*x^6*y^5*z^5*R+2268*x^6*y^5*z^4*R+2268*x^6*y^4*z^5*R-2548*x^6*y^4*z^4*R+2592*x^5*y^5*z^5*R-2916*x^5*y^5*z^4*R-2916*x^5*y^4*z^5*R+864*x^5*y^3*z^3*R-480*x^4*y^3*z^3*R-672*x^6*y^3*z^3*R+192*x^7*y^3*z^3*R-1820*x^4*y^4*z^4*R+336*x^7*y^5*z^3*R-376*x^7*y^4*z^3*R-1176*x^6*y^5*z^3*R+1316*x^6*y^4*z^3*R+1512*x^5*y^5*z^3*R-324*x^5*y^2*z^5*R+168*x^3*y^5*z^3*R-840*x^4*y^5*z^3*R-324*x^3*y^4*z^5*R-324*x^3*y^5*z^4*R-6*x^3*z+3*x^4*z-12*x^3*z^3+6*x^4*z^3+18*x^3*z^2-9*x^4*z^2+6*x^2*y^3-12*x^3*y^3+6*x^4*y^3+18*x^3*y^2-9*x^4*y^2-6*x^3*y+3*x^4*y+288*x^3*y^5*z^5*R+4*x^3*y^2*z^2*R+180*x^4*y^2*z^5*R-864*x^5*y^6*z^5*R+672*x^6*y^6*z^5*R-756*x^6*y^6*z^4*R+280*x^4*y^6*z^3*R+392*x^6*y^6*z^3*R-504*x^5*y^6*z^3*R+480*x^4*y^6*z^5*R-540*x^4*y^6*z^4*R+40*x^3*y^4*z^2*R+40*x^3*y^2*z^4*R-56*x^3*y^6*z^3*R+108*x^3*y^6*z^4*R-96*x^3*y^6*z^5*R+972*x^5*y^6*z^4*R+3276*x^5*y^4*z^4*R+108*x^5*y^2*z^6*R-60*x^4*y^2*z^6*R+140*x^6*y^2*z^3*R-40*x^7*y^2*z^3*R+252*x^6*y^2*z^5*R-72*x^7*y^2*z^5*R-280*x^6*y^2*z^4*R+80*x^7*y^2*z^4*R-324*x^5*y^5*z^2*R+252*x^6*y^5*z^2*R-72*x^7*y^5*z^2*R-280*x^6*y^4*z^2*R+80*x^7*y^4*z^2*R+140*x^6*y^3*z^2*R-40*x^7*y^3*z^2*R-36*x^3*y^5*z^2*R+972*x^5*y^4*z^6*R-864*x^5*y^5*z^6*R-756*x^6*y^4*z^6*R+672*x^6*y^5*z^6*R-504*x^5*y^3*z^6*R+392*x^6*y^3*z^6*R+288*x^5*y^6*z^6*R-224*x^6*y^6*z^6*R+12*x^3*y^6*z^2*R-28*x^6*y^2*z^2*R+32*x^3*y^6*z^6*R+180*x^4*y^5*z^2*R-60*x^4*y^6*z^2*R+108*x^5*y^6*z^2*R-84*x^6*y^6*z^2*R-84*x^6*y^2*z^6*R+8*x^7*y^2*z^2*R-192*x^7*y^6*z^5*R+216*x^7*y^6*z^4*R-112*x^7*y^6*z^3*R+216*x^7*y^4*z^6*R-192*x^7*y^5*z^6*R-112*x^7*y^3*z^6*R+64*x^7*y^6*z^6*R+24*x^7*y^6*z^2*R+24*x^7*y^2*z^6*R-1692*x^5*y^3*z^4*R-188*x^3*y^3*z^4*R-200*x^4*y^2*z^4*R+360*x^5*y^2*z^4*R+364*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+100*x^4*y^3*z^2*R-180*x^5*y^3*z^2*R+96*x^3*y^3*z^3*R+100*x^4*y^2*z^3*R-180*x^5*y^2*z^3*R-200*x^4*y^4*z^2*R+360*x^5*y^4*z^2*R-188*x^3*y^4*z^3*R+12*x^3*y^2*z^6*R-540*x^4*y^4*z^6*R+480*x^4*y^5*z^6*R+280*x^4*y^3*z^6*R+108*x^3*y^4*z^6*R-96*x^3*y^5*z^6*R-56*x^3*y^3*z^6*R+36*x^2*y*z^2+3*y*z^2-9*x^2*z^2-9*y^2*z^2+3*x^2*y-y*z+3*x^2*z-9*x^2*y^2+3*y^2*z)/R
			,
			3*y*y+ -- dy p term
		    (18*x*y*z^2-24*x^3*y*z^3+36*x^3*y*z^2+36*x^2*y*z^3-36*x^2*y^2*z+18*x*y^2*z+24*x^3*y^2*z^3-36*x^3*y^2*z^2-36*x^2*y^2*z^3+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+24*x^3*y^2*z+36*x^2*y^3*z-6*x^2*z^3+6*y^3*z-24*x*y^3*z-6*y^2*z^3-18*y^3*z^2+12*y^3*z^3+24*x*y^2*z^3+36*x*y^3*z^2-24*x*y^3*z^3-12*x*y*z^3+2*x*z^3-1620*x^4*y^5*z^5*R+1638*x^4*y^5*z^4*R+12*x*y^4*z+12*x^3*y^4*z-18*x^2*y^4*z+12*x*y^4*z^3-18*x*y^4*z^2-3*y^4*z+120*x^3*y^3*z^5*R+9*y^4*z^2-6*y^4*z^3-180*x^4*y^3*z^5*R+900*x^4*y^4*z^5*R+54*x^2*y^2*z^2-4*x^3*y^3*z^2*R-420*x^4*y^6*z^6*R+182*x^4*y^3*z^4*R-48*x^6*y^3*z^5*R+48*x^6*y^3*z^4*R+144*x^5*y^3*z^5*R-240*x^5*y^4*z^3*R+320*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+432*x^6*y^5*z^4*R+240*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-1296*x^5*y^5*z^4*R-720*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-144*x^6*y^5*z^3*R+80*x^6*y^4*z^3*R+432*x^5*y^5*z^3*R+432*x^3*y^5*z^3*R-576*x^4*y^5*z^3*R-600*x^3*y^4*z^5*R-1116*x^3*y^5*z^4*R+2*x^3*z+4*x^3*z^3-6*x^3*z^2-18*x^2*y^3+12*x^3*y^3-6*x^3*y^2+1080*x^3*y^5*z^5*R+6*y^3*x-3*y^4*x+9*y^4*x^2-6*y^4*x^3-1008*x^5*y^6*z^5*R+336*x^6*y^6*z^5*R-336*x^6*y^6*z^4*R+448*x^4*y^6*z^3*R+112*x^6*y^6*z^3*R-336*x^5*y^6*z^3*R+1260*x^4*y^6*z^5*R-1274*x^4*y^6*z^4*R+20*x^3*y^4*z^2*R-336*x^3*y^6*z^3*R+868*x^3*y^6*z^4*R-840*x^3*y^6*z^5*R+1008*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-36*x^3*y^5*z^2*R+240*x^5*y^4*z^6*R-432*x^5*y^5*z^6*R-80*x^6*y^4*z^6*R+144*x^6*y^5*z^6*R-48*x^5*y^3*z^6*R+16*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+28*x^3*y^6*z^2*R+280*x^3*y^6*z^6*R+18*x^4*y^5*z^2*R-14*x^4*y^6*z^2*R-144*x^5*y^3*z^4*R-124*x^3*y^3*z^4*R+620*x^3*y^4*z^4*R+2*x^4*y^3*z^2*R+48*x^3*y^3*z^3*R-10*x^4*y^4*z^2*R-240*x^3*y^4*z^3*R-300*x^4*y^4*z^6*R+540*x^4*y^5*z^6*R+60*x^4*y^3*z^6*R+200*x^3*y^4*z^6*R-360*x^3*y^5*z^6*R-40*x^3*y^3*z^6*R-54*x^2*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*z-3*x^2*z-3*x*y^2+9*x^2*y^2-3*y^2*z-10*x^2*y^4*z^2*R-16*x^2*y^3*z^3*R-80*x^3*y^7*z^6*R+2*x^2*y^3*z^2*R+112*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-266*x^2*y^6*z^4*R+252*x^2*y^6*z^5*R+180*x^2*y^4*z^5*R+364*x^4*y^7*z^4*R-360*x^4*y^7*z^5*R+288*x^5*y^7*z^5*R-288*x^5*y^7*z^4*R+96*x^3*y^7*z^3*R+96*x^5*y^7*z^3*R-128*x^4*y^7*z^3*R+240*x^3*y^7*z^5*R-248*x^3*y^7*z^4*R+18*x^2*y^5*z^2*R+38*x^2*y^3*z^4*R-32*x^2*y^7*z^3*R+76*x^2*y^7*z^4*R-72*x^2*y^7*z^5*R-14*x^2*y^6*z^2*R+120*x^4*y^7*z^6*R-96*x^5*y^7*z^6*R+4*x^2*y^7*z^2*R+24*x^2*y^7*z^6*R-8*x^3*y^7*z^2*R+4*x^4*y^7*z^2*R-96*x^6*y^7*z^5*R+96*x^6*y^7*z^4*R-32*x^6*y^7*z^3*R+32*x^6*y^7*z^6*R-190*x^2*y^4*z^4*R+342*x^2*y^5*z^4*R-36*x^2*y^3*z^5*R+80*x^2*y^4*z^3*R-144*x^2*y^5*z^3*R+12*x^2*y^3*z^6*R+108*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-60*x^2*y^4*z^6*R)/R
		   	,
			3*z*z + -- dz p term
			(18*x*y*z^2-24*x^3*y*z^3+24*x^3*y*z^2+36*x^2*y*z^3-54*x^2*y^2*z+18*x*y^2*z+24*x^3*y^3*z^2-36*x^3*y^2*z^2-36*x^2*y^3*z^2+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+36*x^3*y^2*z+36*x^2*y^3*z-18*x^2*z^3-12*x*y^3*z-18*y^2*z^3-6*y^3*z^2+12*y^3*z^3+36*x*y^2*z^3+24*x*y^3*z^2-24*x*y^3*z^3+6*y*z^3-24*x*y*z^3+6*x*z^3-1620*x^4*y^5*z^5*R+900*x^4*y^5*z^4*R+432*x^3*y^3*z^5*R-576*x^4*y^3*z^5*R+1638*x^4*y^4*z^5*R+54*x^2*y^2*z^2+9*x^2*z^4-3*x*z^4-4*x^3*y^2*z^3*R-420*x^4*y^6*z^6*R+320*x^4*y^3*z^4*R-144*x^6*y^3*z^5*R+80*x^6*y^3*z^4*R+432*x^5*y^3*z^5*R-144*x^5*y^4*z^3*R+182*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+240*x^6*y^5*z^4*R+432*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-720*x^5*y^5*z^4*R-1296*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-48*x^6*y^5*z^3*R+48*x^6*y^4*z^3*R+144*x^5*y^5*z^3*R+120*x^3*y^5*z^3*R-180*x^4*y^5*z^3*R-1116*x^3*y^4*z^5*R-600*x^3*y^5*z^4*R+12*x^3*z^3-6*x^3*z^2-6*x^2*y^3+4*x^3*y^3-6*x^3*y^2+2*x^3*y+1080*x^3*y^5*z^5*R-18*x^2*y*z^4+12*x^3*y*z^4+2*y^3*x-6*z^4*x^3+18*x^4*y^2*z^5*R-432*x^5*y^6*z^5*R+144*x^6*y^6*z^5*R-80*x^6*y^6*z^4*R+60*x^4*y^6*z^3*R+16*x^6*y^6*z^3*R-48*x^5*y^6*z^3*R+540*x^4*y^6*z^5*R-300*x^4*y^6*z^4*R+20*x^3*y^2*z^4*R-40*x^3*y^6*z^3*R+200*x^3*y^6*z^4*R-360*x^3*y^6*z^5*R+240*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-14*x^4*y^2*z^6*R+1008*x^5*y^4*z^6*R-1008*x^5*y^5*z^6*R-336*x^6*y^4*z^6*R+336*x^6*y^5*z^6*R-336*x^5*y^3*z^6*R+112*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+280*x^3*y^6*z^6*R-240*x^5*y^3*z^4*R-240*x^3*y^3*z^4*R-10*x^4*y^2*z^4*R+620*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+48*x^3*y^3*z^3*R+2*x^4*y^2*z^3*R-124*x^3*y^4*z^3*R+28*x^3*y^2*z^6*R-1274*x^4*y^4*z^6*R+1260*x^4*y^5*z^6*R+448*x^4*y^3*z^6*R+868*x^3*y^4*z^6*R-840*x^3*y^5*z^6*R-336*x^3*y^3*z^6*R-36*x^2*y*z^2-3*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*y-3*x^2*y-3*x*y^2+9*x^2*y^2-18*x*y^2*z^4+12*x*y^3*z^4+12*x*y*z^4+9*y^2*z^4-6*y^3*z^4-3*y*z^4-16*x^2*y^3*z^3*R+12*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-60*x^2*y^6*z^4*R+108*x^2*y^6*z^5*R+342*x^2*y^4*z^5*R+80*x^2*y^3*z^4*R-190*x^2*y^4*z^4*R+180*x^2*y^5*z^4*R-144*x^2*y^3*z^5*R+38*x^2*y^4*z^3*R-36*x^2*y^5*z^3*R+112*x^2*y^3*z^6*R+252*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-266*x^2*y^4*z^6*R+2*x^2*y^2*z^3*R-10*x^2*y^2*z^4*R-80*x^3*y^6*z^7*R+18*x^2*y^2*z^5*R+4*x^4*y^2*z^7*R-8*x^3*y^2*z^7*R+364*x^4*y^4*z^7*R-360*x^4*y^5*z^7*R-288*x^5*y^4*z^7*R+288*x^5*y^5*z^7*R-128*x^4*y^3*z^7*R+96*x^5*y^3*z^7*R+120*x^4*y^6*z^7*R-96*x^5*y^6*z^7*R+24*x^2*y^6*z^7*R+96*x^6*y^4*z^7*R-96*x^6*y^5*z^7*R-32*x^6*y^3*z^7*R+32*x^6*y^6*z^7*R-14*x^2*y^2*z^6*R+4*x^2*y^2*z^7*R-248*x^3*y^4*z^7*R+240*x^3*y^5*z^7*R+96*x^3*y^3*z^7*R+76*x^2*y^4*z^7*R-72*x^2*y^5*z^7*R-32*x^2*y^3*z^7*R)/R	
		end	
	else
		function source3d(x, y,z, t)
			return 
			2*x +  -- dx p term
		    2*(-18*x*y*z^2-12*x^4*y*z^3+18*x^4*y*z^2+24*x^3*y*z^3-36*x^3*y*z^2-24*x^2*y*z^3+36*x^2*y^2*z-18*x*y^2*z-24*x^2*y^3*z^3+36*x^2*y^3*z^2+36*x^2*y^2*z^3-18*x^2*y*z+6*x*y*z+24*x^3*y*z-12*x^4*y*z+54*x*y^2*z^2-12*x^4*y^3*z+18*x^4*y^2*z+24*x^3*y^3*z-36*x^3*y^2*z-24*x^2*y^3*z+6*x^2*z^3-2*y^3*z+12*x*y^3*z+6*y^2*z^3+6*y^3*z^2-4*y^3*z^3-36*x*y^2*z^3-36*x*y^3*z^2+24*x*y^3*z^3-2*y*z^3+12*x*y*z^3-1440*x^4*y^5*z^5*R+1620*x^4*y^5*z^4*R+168*x^3*y^3*z^5*R-840*x^4*y^3*z^5*R+1620*x^4*y^4*z^5*R-54*x^2*y^2*z^2+36*x^5*y^2*z^2*R-20*x^4*y^2*z^2*R-20*x^3*y^3*z^2*R-20*x^3*y^2*z^3*R-160*x^4*y^6*z^6*R+940*x^4*y^3*z^4*R+336*x^7*y^3*z^5*R-376*x^7*y^3*z^4*R-1176*x^6*y^3*z^5*R+1316*x^6*y^3*z^4*R+1512*x^5*y^3*z^5*R-1692*x^5*y^4*z^3*R+940*x^4*y^4*z^3*R+576*x^7*y^5*z^5*R-648*x^7*y^5*z^4*R-648*x^7*y^4*z^5*R+728*x^7*y^4*z^4*R-2016*x^6*y^5*z^5*R+2268*x^6*y^5*z^4*R+2268*x^6*y^4*z^5*R-2548*x^6*y^4*z^4*R+2592*x^5*y^5*z^5*R-2916*x^5*y^5*z^4*R-2916*x^5*y^4*z^5*R+864*x^5*y^3*z^3*R-480*x^4*y^3*z^3*R-672*x^6*y^3*z^3*R+192*x^7*y^3*z^3*R-1820*x^4*y^4*z^4*R+336*x^7*y^5*z^3*R-376*x^7*y^4*z^3*R-1176*x^6*y^5*z^3*R+1316*x^6*y^4*z^3*R+1512*x^5*y^5*z^3*R-324*x^5*y^2*z^5*R+168*x^3*y^5*z^3*R-840*x^4*y^5*z^3*R-324*x^3*y^4*z^5*R-324*x^3*y^5*z^4*R-6*x^3*z+3*x^4*z-12*x^3*z^3+6*x^4*z^3+18*x^3*z^2-9*x^4*z^2+6*x^2*y^3-12*x^3*y^3+6*x^4*y^3+18*x^3*y^2-9*x^4*y^2-6*x^3*y+3*x^4*y+288*x^3*y^5*z^5*R+4*x^3*y^2*z^2*R+180*x^4*y^2*z^5*R-864*x^5*y^6*z^5*R+672*x^6*y^6*z^5*R-756*x^6*y^6*z^4*R+280*x^4*y^6*z^3*R+392*x^6*y^6*z^3*R-504*x^5*y^6*z^3*R+480*x^4*y^6*z^5*R-540*x^4*y^6*z^4*R+40*x^3*y^4*z^2*R+40*x^3*y^2*z^4*R-56*x^3*y^6*z^3*R+108*x^3*y^6*z^4*R-96*x^3*y^6*z^5*R+972*x^5*y^6*z^4*R+3276*x^5*y^4*z^4*R+108*x^5*y^2*z^6*R-60*x^4*y^2*z^6*R+140*x^6*y^2*z^3*R-40*x^7*y^2*z^3*R+252*x^6*y^2*z^5*R-72*x^7*y^2*z^5*R-280*x^6*y^2*z^4*R+80*x^7*y^2*z^4*R-324*x^5*y^5*z^2*R+252*x^6*y^5*z^2*R-72*x^7*y^5*z^2*R-280*x^6*y^4*z^2*R+80*x^7*y^4*z^2*R+140*x^6*y^3*z^2*R-40*x^7*y^3*z^2*R-36*x^3*y^5*z^2*R+972*x^5*y^4*z^6*R-864*x^5*y^5*z^6*R-756*x^6*y^4*z^6*R+672*x^6*y^5*z^6*R-504*x^5*y^3*z^6*R+392*x^6*y^3*z^6*R+288*x^5*y^6*z^6*R-224*x^6*y^6*z^6*R+12*x^3*y^6*z^2*R-28*x^6*y^2*z^2*R+32*x^3*y^6*z^6*R+180*x^4*y^5*z^2*R-60*x^4*y^6*z^2*R+108*x^5*y^6*z^2*R-84*x^6*y^6*z^2*R-84*x^6*y^2*z^6*R+8*x^7*y^2*z^2*R-192*x^7*y^6*z^5*R+216*x^7*y^6*z^4*R-112*x^7*y^6*z^3*R+216*x^7*y^4*z^6*R-192*x^7*y^5*z^6*R-112*x^7*y^3*z^6*R+64*x^7*y^6*z^6*R+24*x^7*y^6*z^2*R+24*x^7*y^2*z^6*R-1692*x^5*y^3*z^4*R-188*x^3*y^3*z^4*R-200*x^4*y^2*z^4*R+360*x^5*y^2*z^4*R+364*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+100*x^4*y^3*z^2*R-180*x^5*y^3*z^2*R+96*x^3*y^3*z^3*R+100*x^4*y^2*z^3*R-180*x^5*y^2*z^3*R-200*x^4*y^4*z^2*R+360*x^5*y^4*z^2*R-188*x^3*y^4*z^3*R+12*x^3*y^2*z^6*R-540*x^4*y^4*z^6*R+480*x^4*y^5*z^6*R+280*x^4*y^3*z^6*R+108*x^3*y^4*z^6*R-96*x^3*y^5*z^6*R-56*x^3*y^3*z^6*R+36*x^2*y*z^2+3*y*z^2-9*x^2*z^2-9*y^2*z^2+3*x^2*y-y*z+3*x^2*z-9*x^2*y^2+3*y^2*z)/R
			,
			2*y + -- dy p term
		   (18*x*y*z^2-24*x^3*y*z^3+36*x^3*y*z^2+36*x^2*y*z^3-36*x^2*y^2*z+18*x*y^2*z+24*x^3*y^2*z^3-36*x^3*y^2*z^2-36*x^2*y^2*z^3+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+24*x^3*y^2*z+36*x^2*y^3*z-6*x^2*z^3+6*y^3*z-24*x*y^3*z-6*y^2*z^3-18*y^3*z^2+12*y^3*z^3+24*x*y^2*z^3+36*x*y^3*z^2-24*x*y^3*z^3-12*x*y*z^3+2*x*z^3-1620*x^4*y^5*z^5*R+1638*x^4*y^5*z^4*R+12*x*y^4*z+12*x^3*y^4*z-18*x^2*y^4*z+12*x*y^4*z^3-18*x*y^4*z^2-3*y^4*z+120*x^3*y^3*z^5*R+9*y^4*z^2-6*y^4*z^3-180*x^4*y^3*z^5*R+900*x^4*y^4*z^5*R+54*x^2*y^2*z^2-4*x^3*y^3*z^2*R-420*x^4*y^6*z^6*R+182*x^4*y^3*z^4*R-48*x^6*y^3*z^5*R+48*x^6*y^3*z^4*R+144*x^5*y^3*z^5*R-240*x^5*y^4*z^3*R+320*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+432*x^6*y^5*z^4*R+240*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-1296*x^5*y^5*z^4*R-720*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-144*x^6*y^5*z^3*R+80*x^6*y^4*z^3*R+432*x^5*y^5*z^3*R+432*x^3*y^5*z^3*R-576*x^4*y^5*z^3*R-600*x^3*y^4*z^5*R-1116*x^3*y^5*z^4*R+2*x^3*z+4*x^3*z^3-6*x^3*z^2-18*x^2*y^3+12*x^3*y^3-6*x^3*y^2+1080*x^3*y^5*z^5*R+6*y^3*x-3*y^4*x+9*y^4*x^2-6*y^4*x^3-1008*x^5*y^6*z^5*R+336*x^6*y^6*z^5*R-336*x^6*y^6*z^4*R+448*x^4*y^6*z^3*R+112*x^6*y^6*z^3*R-336*x^5*y^6*z^3*R+1260*x^4*y^6*z^5*R-1274*x^4*y^6*z^4*R+20*x^3*y^4*z^2*R-336*x^3*y^6*z^3*R+868*x^3*y^6*z^4*R-840*x^3*y^6*z^5*R+1008*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-36*x^3*y^5*z^2*R+240*x^5*y^4*z^6*R-432*x^5*y^5*z^6*R-80*x^6*y^4*z^6*R+144*x^6*y^5*z^6*R-48*x^5*y^3*z^6*R+16*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+28*x^3*y^6*z^2*R+280*x^3*y^6*z^6*R+18*x^4*y^5*z^2*R-14*x^4*y^6*z^2*R-144*x^5*y^3*z^4*R-124*x^3*y^3*z^4*R+620*x^3*y^4*z^4*R+2*x^4*y^3*z^2*R+48*x^3*y^3*z^3*R-10*x^4*y^4*z^2*R-240*x^3*y^4*z^3*R-300*x^4*y^4*z^6*R+540*x^4*y^5*z^6*R+60*x^4*y^3*z^6*R+200*x^3*y^4*z^6*R-360*x^3*y^5*z^6*R-40*x^3*y^3*z^6*R-54*x^2*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*z-3*x^2*z-3*x*y^2+9*x^2*y^2-3*y^2*z-10*x^2*y^4*z^2*R-16*x^2*y^3*z^3*R-80*x^3*y^7*z^6*R+2*x^2*y^3*z^2*R+112*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-266*x^2*y^6*z^4*R+252*x^2*y^6*z^5*R+180*x^2*y^4*z^5*R+364*x^4*y^7*z^4*R-360*x^4*y^7*z^5*R+288*x^5*y^7*z^5*R-288*x^5*y^7*z^4*R+96*x^3*y^7*z^3*R+96*x^5*y^7*z^3*R-128*x^4*y^7*z^3*R+240*x^3*y^7*z^5*R-248*x^3*y^7*z^4*R+18*x^2*y^5*z^2*R+38*x^2*y^3*z^4*R-32*x^2*y^7*z^3*R+76*x^2*y^7*z^4*R-72*x^2*y^7*z^5*R-14*x^2*y^6*z^2*R+120*x^4*y^7*z^6*R-96*x^5*y^7*z^6*R+4*x^2*y^7*z^2*R+24*x^2*y^7*z^6*R-8*x^3*y^7*z^2*R+4*x^4*y^7*z^2*R-96*x^6*y^7*z^5*R+96*x^6*y^7*z^4*R-32*x^6*y^7*z^3*R+32*x^6*y^7*z^6*R-190*x^2*y^4*z^4*R+342*x^2*y^5*z^4*R-36*x^2*y^3*z^5*R+80*x^2*y^4*z^3*R-144*x^2*y^5*z^3*R+12*x^2*y^3*z^6*R+108*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-60*x^2*y^4*z^6*R)/R
		   	,
		     2*z + -- dz p term
		    (18*x*y*z^2-24*x^3*y*z^3+24*x^3*y*z^2+36*x^2*y*z^3-54*x^2*y^2*z+18*x*y^2*z+24*x^3*y^3*z^2-36*x^3*y^2*z^2-36*x^2*y^3*z^2+18*x^2*y*z-6*x*y*z-12*x^3*y*z-36*x*y^2*z^2-24*x^3*y^3*z+36*x^3*y^2*z+36*x^2*y^3*z-18*x^2*z^3-12*x*y^3*z-18*y^2*z^3-6*y^3*z^2+12*y^3*z^3+36*x*y^2*z^3+24*x*y^3*z^2-24*x*y^3*z^3+6*y*z^3-24*x*y*z^3+6*x*z^3-1620*x^4*y^5*z^5*R+900*x^4*y^5*z^4*R+432*x^3*y^3*z^5*R-576*x^4*y^3*z^5*R+1638*x^4*y^4*z^5*R+54*x^2*y^2*z^2+9*x^2*z^4-3*x*z^4-4*x^3*y^2*z^3*R-420*x^4*y^6*z^6*R+320*x^4*y^3*z^4*R-144*x^6*y^3*z^5*R+80*x^6*y^3*z^4*R+432*x^5*y^3*z^5*R-144*x^5*y^4*z^3*R+182*x^4*y^4*z^3*R-432*x^6*y^5*z^5*R+240*x^6*y^5*z^4*R+432*x^6*y^4*z^5*R-240*x^6*y^4*z^4*R+1296*x^5*y^5*z^5*R-720*x^5*y^5*z^4*R-1296*x^5*y^4*z^5*R+48*x^5*y^3*z^3*R-64*x^4*y^3*z^3*R-16*x^6*y^3*z^3*R-910*x^4*y^4*z^4*R-48*x^6*y^5*z^3*R+48*x^6*y^4*z^3*R+144*x^5*y^5*z^3*R+120*x^3*y^5*z^3*R-180*x^4*y^5*z^3*R-1116*x^3*y^4*z^5*R-600*x^3*y^5*z^4*R+12*x^3*z^3-6*x^3*z^2-6*x^2*y^3+4*x^3*y^3-6*x^3*y^2+2*x^3*y+1080*x^3*y^5*z^5*R-18*x^2*y*z^4+12*x^3*y*z^4+2*y^3*x-6*z^4*x^3+18*x^4*y^2*z^5*R-432*x^5*y^6*z^5*R+144*x^6*y^6*z^5*R-80*x^6*y^6*z^4*R+60*x^4*y^6*z^3*R+16*x^6*y^6*z^3*R-48*x^5*y^6*z^3*R+540*x^4*y^6*z^5*R-300*x^4*y^6*z^4*R+20*x^3*y^2*z^4*R-40*x^3*y^6*z^3*R+200*x^3*y^6*z^4*R-360*x^3*y^6*z^5*R+240*x^5*y^6*z^4*R+720*x^5*y^4*z^4*R-14*x^4*y^2*z^6*R+1008*x^5*y^4*z^6*R-1008*x^5*y^5*z^6*R-336*x^6*y^4*z^6*R+336*x^6*y^5*z^6*R-336*x^5*y^3*z^6*R+112*x^6*y^3*z^6*R+336*x^5*y^6*z^6*R-112*x^6*y^6*z^6*R+280*x^3*y^6*z^6*R-240*x^5*y^3*z^4*R-240*x^3*y^3*z^4*R-10*x^4*y^2*z^4*R+620*x^3*y^4*z^4*R-36*x^3*y^2*z^5*R+48*x^3*y^3*z^3*R+2*x^4*y^2*z^3*R-124*x^3*y^4*z^3*R+28*x^3*y^2*z^6*R-1274*x^4*y^4*z^6*R+1260*x^4*y^5*z^6*R+448*x^4*y^3*z^6*R+868*x^3*y^4*z^6*R-840*x^3*y^5*z^6*R-336*x^3*y^3*z^6*R-36*x^2*y*z^2-3*y*z^2-3*x*z^2+9*x^2*z^2+9*y^2*z^2+x*y-3*x^2*y-3*x*y^2+9*x^2*y^2-18*x*y^2*z^4+12*x*y^3*z^4+12*x*y*z^4+9*y^2*z^4-6*y^3*z^4-3*y*z^4-16*x^2*y^3*z^3*R+12*x^2*y^6*z^3*R-324*x^2*y^5*z^5*R-60*x^2*y^6*z^4*R+108*x^2*y^6*z^5*R+342*x^2*y^4*z^5*R+80*x^2*y^3*z^4*R-190*x^2*y^4*z^4*R+180*x^2*y^5*z^4*R-144*x^2*y^3*z^5*R+38*x^2*y^4*z^3*R-36*x^2*y^5*z^3*R+112*x^2*y^3*z^6*R+252*x^2*y^5*z^6*R-84*x^2*y^6*z^6*R-266*x^2*y^4*z^6*R+2*x^2*y^2*z^3*R-10*x^2*y^2*z^4*R-80*x^3*y^6*z^7*R+18*x^2*y^2*z^5*R+4*x^4*y^2*z^7*R-8*x^3*y^2*z^7*R+364*x^4*y^4*z^7*R-360*x^4*y^5*z^7*R-288*x^5*y^4*z^7*R+288*x^5*y^5*z^7*R-128*x^4*y^3*z^7*R+96*x^5*y^3*z^7*R+120*x^4*y^6*z^7*R-96*x^5*y^6*z^7*R+24*x^2*y^6*z^7*R+96*x^6*y^4*z^7*R-96*x^6*y^5*z^7*R-32*x^6*y^3*z^7*R+32*x^6*y^6*z^7*R-14*x^2*y^2*z^6*R+4*x^2*y^2*z^7*R-248*x^3*y^4*z^7*R+240*x^3*y^5*z^7*R+96*x^3*y^3*z^7*R+76*x^2*y^4*z^7*R-72*x^2*y^5*z^7*R-32*x^2*y^3*z^7*R)/R	
		end	
	end
	
end

function inletVel2d(x, y, t)
	return usol2d(x, y, t), vsol2d(x, y, t)
end

function inletVel3d(x, y, z, t)
	return usol3d(x,y,z,t), vsol3d(x,y,z,t), wsol3d(x,y,z,t)
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

function CreateDomainDisc(approxSpace, discType, p)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity(1.0/R);
	NavierStokesDisc:set_source("source"..dim.."d")
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(bPac)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		NavierStokesDisc:set_stabilization(3)
	end
	
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Boundary")
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	
	requiredSubsets = {"Inner", "Boundary"}
	local dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
	
	return dom
end

function CreateApproxSpace(dom, discType, p)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, p, p-1)
	
	-- print statistic on the distributed dofs
	--approxSpace:init_levels()
	--approxSpace:init_top_surface()
	--approxSpace:print_statistic()
	--approxSpace:print_local_dof_statistic(2)
	
	return approxSpace
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

function CreateSolver(approxSpace, discType, p)

	local base = nil
	if discType == "fvcr" then
		base =  LinearSolver()
		base:set_preconditioner(DiagVanka())
		base:set_convergence_check(ConvCheck(10000, 1e-7, 1e-3, false))
	else
		base = SuperLU()
	end
	
	local smoother = nil
	if discType == "fvcr" then 
		smoother = Vanka()
	else
		local smooth = util.smooth.parseParams()
		smoother = util.smooth.create(smooth)
	end
	
	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	--gmg:set_damp(MinimalResiduumDamping())
	--gmg:set_damp(MinimalEnergyDamping())
	gmg:add_prolongation_post_process(AverageComponent("p"))
	--gmg:set_debug(dbgWriter)
	
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	solver:set_convergence_check(ConvCheck(100000, 1e-10, 1e-99, false))
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(ConvCheck(100, nlintol, nlinred, false))
	--newtonSolver:set_line_search(StandardLineSearch(30, 1.0, 0.85, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)

	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
	AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

if not(bInstat) then

	if bConvRates then
		
		local options = {	
		
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
		}
	
		if util.HasParamOption("-replot") then
			util.rates.static.replot(options); exit()
		end
	
		util.rates.static.compute(
		{
			ExactSol = {
				["u"] = "usol"..dim.."d",
				["v"] = "vsol"..dim.."d",
				["p"] = "psol"..dim.."d"
			},
			ExactGrad =  {
				["u"] = "ugrad"..dim.."d",
				["v"] = "vgrad"..dim.."d",
				["p"] = "pgrad"..dim.."d"
			},
			
			PlotCmps = { v = {"u","v"}, p = {"p"}},
			MeasLabel = function (disc, p) return disc.." $\\mathbb{Q}_{"..p.."}/\\mathbb{Q}_{"..(p-1).."}$" end,
			
			CreateDomain = CreateDomain,
			CreateApproxSpace = CreateApproxSpace,
			CreateDomainDisc = CreateDomainDisc,
			CreateSolver = CreateSolver,
			
			ComputeSolution = ComputeNonLinearSolution,
			
			DiscTypes = 
			{
			  {type = "fv", pmin = 2, pmax = 4, lmin = 2, lmax = numRefs}
			},
			
			gpOptions = options,
			noplot = true,
		})
	end
	
	if not(bConvRates) then
	
		local p = vorder
		local dom = CreateDomain()
		local approxSpace = CreateApproxSpace(dom, discType, p)
		local domainDisc = CreateDomainDisc(approxSpace, discType, p)
		local solver = CreateSolver(approxSpace, discType, p)
		print(solver:config_string())
		
		local u = GridFunction(approxSpace)
		u:set(0)
		
		timeStart = os.clock()
		ComputeNonLinearSolution(u, domainDisc, solver)
		timeEnd = os.clock()
		print("Computation took " .. timeEnd-timeStart .. " seconds.")
		
		-- to make error computation for p reasonable
		-- p would have to be adjusted by adding a reasonable constant
		local FctCmp = approxSpace:names()
		for d = 1,#FctCmp do
			print("L2Error in '"..FctCmp[d].. "' is ".. 
					L2Error(FctCmp[d].."sol"..dim.."d", u, FctCmp[d], 0.0, 1, "Inner"))
		end
		for d = 1,#FctCmp do
			print("Maximum error in '"..FctCmp[d].. "' is ".. 
					 MaxError(FctCmp[d].."sol"..dim.."d", u, FctCmp[d]))
		end
		
		local VelCmp = {}
		for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
		vtkWriter = VTKOutput()
		vtkWriter:select(VelCmp, "velocity")
		vtkWriter:select("p", "pressure")
		vtkWriter:print("Dirichlet", u)
	end
end

if bInstat then
	step=0;
	time=0;
	
	-- create new grid function for old value
	u = GridFunction(approxSpace)
	u:set(0)
	uOld = u:clone()
	
	vtkWriter = VTKOutput()
	vtkWriter:select(VelCmp, "velocity")
	vtkWriter:select("p", "pressure")
	vtkWriter:print("TimeDirichlet", u, 0,0)
	
	tBefore = os.clock()
	
	-- store grid function in vector of  old solutions
	solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(uOld, time)
	
	for step = 1, numTimeSteps do
		print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

		-- choose time step
		do_dt = dt
	
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
	
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
	
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 

		-- update new time
		time = solTimeSeries:time(0) + do_dt
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()

		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
	
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		-- compute CFL number 
		cflNumber(u,do_dt)

		print("++++++ TIMESTEP " .. step .. "  END ++++++");
		write("\n")
		l2error = L2Error("usol"..dim.."d", u, "u", time, 1, "Inner")
		write("L2Error in u component is "..l2error .."\n")
		l2error = L2Error("vsol"..dim.."d", u, "v", time, 1, "Inner")
		write("L2Error in v component is "..l2error .."\n")
		maxerror = MaxError("usol"..dim.."d", u, "u",time)
		write("Maximum error in u component is "..maxerror .."\n")
		maxerror = MaxError("vsol"..dim.."d", u, "v",time)
		write("Maximum error in v component is "..maxerror .."\n")
		write("\n")
		
		vtkWriter:print("TimeDirichlet", u,step,time)
	end
	
	tAfter = os.clock()
	print("Computation took " .. tAfter-tBefore .. " seconds.");

end