function CreateSolver(approxSpace, discType, p)

	local base = SuperLU()
	
	local smoother = nil
	if discType == "fvcr" or discType == "fecr" then 
		smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	elseif discType == "fv1" then 
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
	end
	
	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	gmg:add_prolongation_post_process(AverageComponent("p"))
	transfer = StdTransfer()
	transfer:enable_p1_lagrange_optimization(false)
	--gmg:set_transfer(transfer)
	
	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, smoother)
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))	
	end
		
	local convCheck = ConvCheck(500, 1e-11, 1e-99, true)
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)

	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
	AdjustMeanValue(u, "p")
end
