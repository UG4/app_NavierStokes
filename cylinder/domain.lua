--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	local dom = Domain()
	LoadDomain(dom, gridName)
	
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner =  GlobalDomainRefiner(dom)
	if     dim == 2 then 
		ProjectVerticesToSphere(dom, {0.2, 0.2}, 0.05, 0.001)
		falloffProjector = SphericalFalloffProjector(dom, {0.2, 0.2}, 0.1, 0.15)
	elseif dim == 3 then 
		falloffProjector = CylindricalFalloffProjector(dom, {0.5, 0.2, 0.0}, {0, 0, 1}, 0.04, 0.1)
	end
	local refProjector = DomainRefinementProjectionHandler(dom)
	refProjector:set_callback("Inner", falloffProjector)
	refProjector:set_callback("CylinderWall", falloffProjector)
	if dim == 3 then
	refProjector:set_callback("BackWall", falloffProjector)
	refProjector:set_callback("FrontWall", falloffProjector)
	end
	refiner:set_refinement_callback(refProjector)
	
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

function CreateApproxSpace(dom, discType, p)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, p, p-1)
	
	-- print statistic on the distributed dofs
	--approxSpace:init_levels()
	--approxSpace:init_top_surface()
	--approxSpace:print_statistic()
	--approxSpace:print_local_dof_statistic(2)
	
	return approxSpace
end