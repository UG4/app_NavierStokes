
-- Create util namespace
util = util or {}


util.ns = util.ns or {}

function util.ns.parseParams()

	local discType = util.GetParam("-type", "fvcr", "Disc type", {"fv1","fv","fe","fecr","fvcr"})
	
	local order, vorder, porder
	if discType == "fe" or discType == "fv" then
		order 		= util.GetParamNumber("-order", 2, "Order of velocity-space, pressure-order accordingly")
		vorder 		= util.GetParamNumber("-vorder", order, "Order of velocity-space")
		porder 		= util.GetParamNumber("-porder", vorder-1, "Order of pressure-space")
	end
	
	return discType, vorder, porder
end

function util.ns.CreateApproxSpace(dom, type, vorder, porder)

	-- create Approx Space w.r.t. to domain
	local dim = dom:get_dim()
	local approxSpace = ApproximationSpace(dom)
	
	-- components of the velocity 
	local VelCmp, FctCmp
	if 		dim == 1 then VelCmp = {"u"}; 			FctCmp = {"u", "p"}
	elseif  dim == 2 then VelCmp = {"u", "v"};	 	FctCmp = {"u", "v", "p"}
	elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"}
	else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end
	
	-- distinguish types
	if type == "fv1" then
		approxSpace:add_fct(VelCmp	, "Lagrange", 1) 
		approxSpace:add_fct("p"		, "Lagrange", 1) 
	elseif type == "fecr" or type == "fvcr" then
		approxSpace:add_fct(VelCmp	, "Crouzeix-Raviart",1)
		approxSpace:add_fct("p"		, "Piecewise-Constant") 
	elseif type == "fv" or type == "fe" then
		if porder==0 then
			print("ns.createApproximationSpace: porder 0 not valid for "..type); exit();
		end
		approxSpace:add_fct(VelCmp	, "Lagrange", vorder) 
		approxSpace:add_fct("p"		, "Lagrange", porder) 
	else 
		print("ns.createApproximationSpace: Disc Type '"..type.."' not supported."); exit(); 
	end
	
	return approxSpace, FctCmp, VelCmp
end


util.gmg = util.gmg or {}

function util.gmg.parseParams()

	local numPreSmooth, numPostSmooth
	if util.HasParamOption("-numSmooth") then
		numPreSmooth  = util.GetParamNumber("-numSmooth", 2, "Number pre/post-smoothing")
		numPostSmooth = numPreSmooth
	else
		numPreSmooth  = util.GetParamNumber("-numPreSmooth", 2, "Number pre-smoothing")
		numPostSmooth = util.GetParamNumber("-numPostSmooth", 2, "Number post-smoothing")
	end
	
	local baseLev = util.GetParamNumber("-baseLev", 0, "Base level")
	local cycle =  util.GetParam("-cycle", "V", "gmg-cycle type", {"V","W","F"})
	local bRAP = util.HasParamOption("-rap", "use rap product as level matrices")
	
	return numPreSmooth, numPostSmooth, baseLev, cycle, bRAP
end

function util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
						 cycle, baseSolver, baseLev, bRAP)

	local gmg = GeometricMultiGrid(approxSpace)
	
	gmg:set_base_level(baseLev)
	gmg:set_base_solver(baseSolver)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	gmg:set_smoother(smoother)
	gmg:set_cycle_type(cycle)
	gmg:set_num_presmooth(numPreSmooth)
	gmg:set_num_postsmooth(numPostSmooth)
	gmg:set_rap(bRAP)
	
	return gmg
end


util.smooth = util.smooth or {}

function util.smooth.parseParams()
	local smooth = util.GetParam("-smooth", "jac", "Smoother Type", 
					{"jac","ilu","ilut","egs","gs","sgs"})
	
	return smooth
end

function util.smooth.create(smooth)

	local smoother = nil
	
	if 	    smooth == "ilu"  then smoother = ILU();
	elseif 	smooth == "ilut" then smoother = ILUT(1e-6);
	elseif 	smooth == "egs"  then smoother = ElementGaussSeidel(groupType);
	elseif 	smooth == "cgs"  then smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
	elseif 	smooth == "jac"  then smoother = Jacobi(0.66);
	elseif 	smooth == "gs"   then smoother = GaussSeidel();
	elseif 	smooth == "sgs"  then smoother = SymmetricGaussSeidel();
	else print("Smoother type '"..smooth.."' not found"); exit(); end
	
	return smoother
end


util.solver = util.solver or {}

function util.solver.parseParams()
	local solver = util.GetParam("-solver", "bicgstab", "Linear Solver Type", 
								{"ls","bicgstab","cg","lu","schur"})
	return solver
end

function util.solver.create(sol, precond)

	local solver = nil
	
	if 		sol == "ls" 		then 
		solver = LinearSolver();
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif 	sol == "bicgstab" 	then 
		solver = BiCGStab();
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif 	sol == "cg" 		then 
		solver = CG();
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif  sol == "schur" then 

		local skeletonSolver = BiCGStab()
		skeletonSolver:set_preconditioner(ILU())
		skeletonSolver:set_convergence_check(ConvCheck(10000, 1e-12, 1e-10, true))	

		local schur = SchurComplement()
		local lu = SuperLU()
		schur:set_dirichlet_solver(lu)
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
	
		solver = LinearSolver()
		solver:set_preconditioner(schur)
			
	elseif sol == "lu" then
		solver = SuperLU()
	else
		print("Solver not found."); exit();
	end	
	
	return solver
end