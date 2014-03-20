--------------------------------------------------------------------------------
--[[!
-- \file apps/navier_stokes/chorin.lua
-- \author Andreas Vogel 
-- \brief Lua - Script to test the time-dependent navier-stokes
]]--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_kinetic.lua")

dim 		= util.GetParamNumber("-dim", 2)
numRefs 	= util.GetParamNumber("-numRefs",4)
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)

startTime  = util.GetParamNumber("-start", 0.0, "start time")
endTime    = util.GetParamNumber("-end", 1e0, "start time")
dt         = util.GetParamNumber("-dt", 1e-1, "time step size")

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
--	gridName = util.GetParam("-grid", "grids/unit_square_01_tri_unstruct_fine.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_1x1.ugx")
else
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_tets.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_1x1x1.ugx")
end
if dim~=2 and dim~=3 then
   print("Chosen Dimension " .. dim .. " not supported. Exiting.") exit() 
end

--------------------------------------------------------------------------------
--  Setup FV Element Discretization
--------------------------------------------------------------------------------

function CreateDomain()
	-- choose algebra
	InitUG(dim, AlgebraType("CPU", 1));
		
	-- Create, Load, Refine and Distribute Domain
	local neededSubsets = {}
	local dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

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
-- Problem
--------------------------------------------------------------------------------
--[[ This problem is adapted after Chorin, 1967
The analytical solution can be constructed e.g. via maple:
with(codegen,C):
# choose divergence free velocity (using g) and pressure p 
g:=-1/s*cos(s*x)*cos(s*y)*exp(-2*s*s*t/R);
u:=-diff(g,y);
v:=diff(g,x); 
p:=-1/4*(cos(2*s*x)+cos(2*s*y))*exp(-4*s*s*t/R);

# rhs is chosen so that (Navier)-Stokes system is fulfilled

time_u := factor(simplify(diff(u,t)));
laplace_u := factor(simplify(-1/R*(diff(u,x,x)+diff(u,y,y)))); 
nonlin_u := factor(simplify(u*diff(u,x)+v*diff(u,y))); 
press_u := factor(simplify(diff(p,x)));

time_v := factor(simplify(diff(v,t)));
laplace_v := factor(simplify(-1/R*(diff(v,x,x)+diff(v,y,y))));
nonlin_v := factor(simplify(u*diff(v,x)+v*diff(v,y)));
press_v := factor(simplify(diff(p,y)));

rhs_u := factor(simplify(time_u + laplace_u + nonlin_u + press_u));
rhs_v := factor(simplify(time_v + laplace_v + nonlin_v + press_v));
C( u ); 
C( v ); 
C( p ); 
C( time_u + laplace_u + nonlin_u + press_u );
C( time_v + laplace_v + nonlin_v + press_v );
--]]

local n = 4
local tau = 1000
local s = n * math.pi

function uSol2d(x, y, t) return -math.cos(s*x)*math.sin(s*y)*math.exp(-2*s*s*t/tau)  end
function vSol2d(x, y, t) return  math.sin(s*x)*math.cos(s*y)*math.exp(-2*s*s*t/tau) end
function pSol2d(x, y, t) return -1/4*(math.cos(2*s*x)+math.cos(2*s*y))*math.exp(-4*s*s*t/tau)  end
--function pSol2d(x, y, t) return 0  end

function uGrad2d(x, y, t) return s*math.sin(s*x)*math.sin(s*y)*math.exp(-2*s*s*t/tau),
								 -s*math.cos(s*x)*math.cos(s*y)*math.exp(-2*s*s*t/tau) end
function vGrad2d(x, y, t) return s*math.cos(s*x)*math.cos(s*y)*math.exp(-2*s*s*t/tau),
								 -s*math.sin(s*x)*math.sin(s*y)*math.exp(-2*s*s*t/tau) end
function pGrad2d(x, y, t) return 2*s/4*(math.sin(2*s*x))*math.exp(-4*s*s*t/tau),
								2*s/4*(math.sin(2*s*y))*math.exp(-4*s*s*t/tau)  end
--function pGrad2d(x, y, t) return 0, 0 end


function inletVel2d(x, y, t)
	return uSol2d(x, y, t), vSol2d(x, y, t)
end

function CreateDomainDisc(approxSpace, discType, p)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity(1.0/tau);
	
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
		--NavierStokesDisc:set_stabilization(3)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(p*p+10)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(p*p+10)
	end
	
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Boundary")
		
	DirichletDisc = DirichletBoundary()
	DirichletDisc:add("uSol2d", "u", "Boundary")
	DirichletDisc:add("vSol2d", "v", "Boundary")
	DirichletDisc:add("pSol2d", "p", "Boundary")
		
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
--	domainDisc:add(DirichletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

function CreateSolver(approxSpace, discType, p)

	local base = nil
	if discType == "fvcr" then
		base =  LinearSolver()
		base:set_preconditioner(DiagVanka())
		base:set_convergence_check(ConvCheck(10000, 1e-7, 1e-3, false))
	else
		base = SuperLU()
		--base = BiCGStab()
		--base:set_preconditioner(ILUT(1e-2))
		--base:set_convergence_check(ConvCheck(10000, 5e-15, 1e-2, true))	
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
	--gmg:set_gathered_base_solver_if_ambiguous(true)
	--gmg:set_rap(true)
	
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-13, 1e-3, true))	
	end
	--solver = SuperLU()
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(ConvCheck(500, 5e-12, 1e-99, true))
	--newtonSolver:set_line_search(StandardLineSearch(30, 1.0, 0.9, true, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

--------------------------------------------------------------------------------
-- Run
--------------------------------------------------------------------------------

gpOpt = 
{
	size = 				{12.5, 8.75}, -- the size of canvas (i.e. plot)
	grid = 				"lc rgb 'grey70' lt 0 lw 1", 

	fontsize =			12,
	fontscale = 		0.8,
	
	datastyle = 		"lines",
	pm3d = 				true,

	grid = 				"lc rgb 'grey70' lt 0 lw 1", 
	linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(5, 1, 360+40, 85, 0.4, 0.4), 
						linewidth = 1, pointsize = 1},
	border = 			" 895 back lc rgb 'grey40' lw 2",
	decimalsign = 		",",

	tics =	 			{x = "border in scale 1 0.01, 2 format '%g \\hspace{0.3cm}' font ',8'",
						 y = "border in scale 1 0.000625, 2 format '\\hspace{1cm} %g' font ',8'", 
						 z = "border in scale 1 format '%.te%01T' font ',8'"}, 

	mtics =	 			true,

	padrange = 			{ x = {0.1,0.1}, y = {0.1,0.1}, z = {0.1,0.1}},
	labeloffset = 		{ x = "offset graph 0,-0.1,0",
						  y = "offset graph 0.1,0,0",
						  z = "offset graph 0.1,0,0.8" },
	key =	 			"on box lc rgb 'grey40' bmargin horizontal maxrows 1 Left reverse spacing 1 width 1.1 samplen 1 height 0.5",
	key = "off",
	"set xrange [] reverse"
}

if util.HasParamOption("-replot") then
	util.rates.kinetic.replot(gpOpt)
	exit()
end

util.rates.kinetic.compute( 
{
	ExactSol = {
		["u"] = "uSol"..dim.."d",
		["v"] = "vSol"..dim.."d",
		["p"] = "pSol"..dim.."d"
	},
	ExactGrad =  {
		["u"] = "uGrad"..dim.."d",
		["v"] = "vGrad"..dim.."d",
		["p"] = "pGrad"..dim.."d"
	},

	PlotCmps = { v = {"u","v"}, p = {"p"}},
	MeasLabel = function (disc, p) return disc.." $\\mathbb{Q}_{"..p.."}/\\mathbb{Q}_{"..(p-1).."}$" end,

	
	CreateDomain = CreateDomain,
	CreateApproxSpace = CreateApproxSpace,
	CreateDomainDisc = CreateDomainDisc,
	CreateSolver = CreateSolver,
	
	StartTime = startTime,
	EndTime = endTime,

	SetStartSolution = function (u, time)
		Interpolate("uSol"..dim.."d", u, "u", time);
		Interpolate("vSol"..dim.."d", u, "v", time);
		Interpolate("pSol"..dim.."d", u, "p", time);
		--AdjustMeanValue(u, "p")
	end,
	
	SpaceDiscs = 
	{
	  {type = "fv", pmin = 2, pmax = 3, lmin = 4, lmax = numRefs} 
	},
	
	
	TimeDiscs =
	{
	  {type = "alexander", orderOrTheta = 3, dt = dt, sub = 2, refs = 0}
	},
	
	gpOptions = gpOpt,
	MaxLevelPadding = util.rates.kinetic.NoMaxLevelPadding,
	best = false,
	noplot = true,
	plotSol = true,	
})
				  