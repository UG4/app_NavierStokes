if not(bConvRates) and not(bBenchmarkRates) then

	local p = vorder
	local dom = CreateDomain()
	local approxSpace = CreateApproxSpace(dom, discType, p)
	local domainDisc = CreateDomainDisc(approxSpace, discType, p)
	local solver = CreateSolver(approxSpace, discType, p)
	--solver:set_debug(GridFunctionDebugWriter(approxSpace))
			
	print(solver:config_string())
	
	local u = GridFunction(approxSpace)
	u:set(0)
	
--	ComputeNonLinearSolution(u, CreateDomainDisc(approxSpace, "fe", p), solver)
	ComputeNonLinearSolution(u, domainDisc, solver)

	local FctCmp = approxSpace:names()
	local VelCmp = {}
	
	for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
	vtkWriter = VTKOutput()
	vtkWriter:select(VelCmp, "velocity")
	vtkWriter:select("p", "pressure")
	vtkWriter:print("Cylinder", u)

	if dim == 2 then
		local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
		local C_D = 2*DL[1]/(Umean2*L)
		local C_L = 2*DL[2]/(Umean2*L)
	
		local PEval = GlobalGridFunctionNumberData(u, "p")
		local Delta_P = PEval:evaluate({0.15, 0.2}) - PEval:evaluate( {0.25, 0.2} )
	
		print("p1: "..PEval:evaluate({0.15, 0.2}))
		print("p2: "..PEval:evaluate({0.25, 0.2}))
	
		print("C_D - ref.CD: "..string.format("%.3e", C_D - ref.CD))
		print("C_L - ref.CL: "..string.format("%.3e", C_L - ref.CL))
		print("Delta_P - ref.DeltaP: "..string.format("%.3e", Delta_P - ref.DeltaP))
	end	
end
