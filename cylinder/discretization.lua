globalNSDisc = nil

function CreateDomainDisc(approxSpace, discType, vorder, porder)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	globalNSDisc = NavierStokesDisc
				
	local porder = approxSpace:lfeid(dim):order()
	local vorder = approxSpace:lfeid(0):order()
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(true)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		NavierStokesDisc:set_stabilization(3)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	
	-- setup Outlet
	--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	--OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		return 4 * Um * y * (H-y) / (H*H), 0.0
	end
	function inletVel3d(x, y, z, t)
		return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0
	end
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
	--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	if dim == 2 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall")
	elseif dim == 3 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall")	
	end
	
	-- Finally we create the discretization object which combines all the
	-- separate discretizations into one domain discretization.
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
	domainDisc:add(WallDisc)
	--domainDisc:add(OutletDisc)
	
	return domainDisc
end
