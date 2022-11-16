
util.transforming = util.transforming or {}
util.transforming.CreateTransformingSmoother = function(approxSpace, VelCmp, type)

 local type = type or "FE"

	local rightTrafoDisc = DomainDiscretization(approxSpace)
	
	for i = 1, #VelCmp do
		print (VelCmp[i])
		local vIdent = DirichletBoundary()
		vIdent:add(ConstUserNumber(0), VelCmp[i], "Inner,Boundary,PressureNode")
		rightTrafoDisc:add(vIdent)
	end
	
	
	local pLaplace = ConvectionDiffusionFE("p", "Inner")
	pLaplace:set_diffusion(1.0)
	pLaplace:set_mass_scale(0.0)
	rightTrafoDisc:add(pLaplace)
	
	
	-- C) Boundary 
	-- TODO: ???
	
	
	-- D) Solver
	local ilu = ILU()
	--ilu:set_sort_type(0)   -- sort: u1, ..., un, v1, ..., vn, p1, .., pn
	--ilu:set_debug(dbgWriter)
	
	--local ilut = ILUT()
	--ilu:set_sort_type(0)   -- sort: u1, ..., un, v1, ..., vn, p1, .., pn
	--ilut:set_debug(dbgWriter)
	
	
	local sgs = SymmetricGaussSeidel()
	sgs:set_damp(1.0)
	
	local jac = Jacobi(0.6)
	
	trafoSmoother = AssembledTransformingSmoother2dCPU1(sgs, rightTrafoDisc)
	
	return trafoSmoother
end 

--[[ 
Creates a Wittum-type transforming smoother:

A) transformedDisc:

-\vec \delta \vec u  = 0
\nabla \cdot u - \beta p =0

where $beta Id\approx  div [-(\vec \delta)^{-1} grad p]$.


B) rightTransformDisc:
   u + \nabla p = 0
             p = p     (pIdent)
--]]

function CreateTransformingSmoother1(approxSpace, VelCmp, beta)

	local rightTrafoDisc = DomainDiscretization(approxSpace)
	local divLinker = ScaleAddLinkerVectorVector2d()

	local unitVecs = {ConstUserVector(0.0), ConstUserVector(0.0)}
	unitVecs[1]:set_entry(1.0, 0.0)
	unitVecs[2]:set_entry(0.0, 1.0)

	-- A) transformedDisc
	local transformedDisc =  DomainDiscretization(approxSpace);
	for i = 1, #VelCmp do
		print (VelCmp[i])
		local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
		tmp:set_diffusion(1.0)
	
		tmp:set_reaction_rate(0.0)
		divLinker:add(unitVecs[i], tmp:gradient())
		transformedDisc:add(tmp)
	end
	local div = ConvectionDiffusion("p", "Inner", type)
	div:set_diffusion(0.0)
	--div:set_vector_source()
	div:set_source(divLinker)
	div:set_reaction_rate(beta)
	transformedDisc:add(div)
	transformedDisc:add(BndDisc)

	-- B) rightTransformDisc
	local rightTrafoDisc = DomainDiscretization(approxSpace)
	for i = 1, #VelCmp do
		local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
		tmp:set_diffusion(0)
		tmp:set_reaction_rate(1)
		tmp:set_reaction(GridFunctionGradientComponentData(u, "p", i))
		rightTrafoDisc:add(tmp)
	end
	pIdent = DirichletBoundary()
	pIdent:add(ConstUserNumber(0), "p", "Inner,Boundary")
	rightTrafoDisc:add(pIdent)
	
	-- C) Boundary 
	-- TODO: ???
	
	
	-- D) Solver
	local ilu = ILU()
	ilu:set_sort_type(3)
	ilu:set_debug(dbgWriter)
	trafoSmoother = AssembledTransformingSmoother2dCPU1(transformedDisc, ilu, rightTrafoDisc)

	return trafoSmoother
end 


--[[
-- This is (supposed to be) Brandt Dinar DGS
-- (Boundary conditions are not 100% clear...)

rightTrafoDisc = DomainDiscretization(approxSpace)
pLaplace = ConvectionDiffusion("p", "Inner", type)
pLaplace:set_diffusion(1)
rightTrafoDisc:add(pLaplace)
OutletDisc = DirichletBoundary()
OutletDisc:add(0.0, "p", "Inlet, Outlet, UpperWall, LowerWall")
rightTrafoDisc:add(OutletDisc)
for i = 1, #VelCmp do
	local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
	tmp:set_diffusion(0)
	tmp:set_reaction_rate(1)
	tmp:set_reaction(GridFunctionGradientComponentData(u, "p", i))
	rightTrafoDisc:add(tmp)
end
tmp = DirichletBoundary()
tmp:add(0.0, "u", "Inlet, Outlet, UpperWall, LowerWall")
tmp:add(0.0, "v", "Inlet, Outlet, UpperWall, LowerWall")
rightTrafoDisc:add(tmp)


TrafoSystenDisc = DomainDiscretization(approxSpace)
for i = 1, #VelCmp do
	local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
	tmp:set_diffusion(Viscosity)
	TrafoSystenDisc:add(tmp)
end
tmp = ConvectionDiffusion("p", "Inner", type)
tmp:set_diffusion(1)
TrafoSystenDisc:add(tmp)
tmp = DirichletBoundary()
tmp:add(0.0, "u", "Inlet, Outlet, UpperWall, LowerWall")
tmp:add(0.0, "v", "Inlet, Outlet, UpperWall, LowerWall")
TrafoSystenDisc:add(tmp)

trafoSmoother = AssembledTransformingSmoother(rightTrafoDisc, TrafoSystemDisc, Jacobi())
trafoSmoother:set_debug(GridFunctionDebugWriter(approxSpace))
trafoSmoother:set_damp(1)
--]]

--[[ 
Creates a Wittum-type transforming smoother:

A) transformedDisc (solve!):

-\vec \delta \vec u          = 0
- \nabla \cdot u     + Q_l p = 0

where $beta Id\approx  div [-(\vec \delta)^{-1} grad p]$.


B) rightTransformDisc (apply!):   
       u:= u + \nabla p 
       p:=      Q^p_l p  
       
--]]

function CreateTransformingSmoother2(approxSpace, VelCmp, beta)

	local rightTrafoDisc = DomainDiscretization(approxSpace)
	local divLinker = ScaleAddLinkerVectorVector2d()

	local unitVecs = {ConstUserVector(0.0), ConstUserVector(0.0)}
	unitVecs[1]:set_entry(1.0, 0.0)
	unitVecs[2]:set_entry(0.0, 1.0)

	-- A) transformedDisc
	local transformedDisc =  DomainDiscretization(approxSpace);
	for i = 1, #VelCmp do
		print (VelCmp[i])
		local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
		tmp:set_diffusion(1)
		tmp:set_reaction_rate(0.0)
		divLinker:add(unitVecs[i], tmp:gradient())
		transformedDisc:add(tmp)
	end
	local div = ConvectionDiffusion("p", "Inner", type)
	div:set_diffusion(1.0)
	div:set_source(divLinker)
	transformedDisc:add(div)
	transformedDisc:add(BndDisc)

	-- B) rightTransformDisc
	local rightTrafoDisc = DomainDiscretization(approxSpace)
	for i = 1, #VelCmp do
		local tmp = ConvectionDiffusion(VelCmp[i], "Inner", type)
		tmp:set_diffusion(0)
		tmp:set_reaction_rate(1)
		tmp:set_reaction(GridFunctionGradientComponentData(u, "p", i))
		rightTrafoDisc:add(tmp)
	end
	pIdent = DirichletBoundary()
	pIdent:add(ConstUserNumber(0), "p", "Inner,Boundary")
	rightTrafoDisc:add(pIdent)
	
	-- C) Boundary 
	-- TODO: ???
	
	
	-- D) Solver
	local ilu = ILU()
	ilu:set_sort_type(3)
	ilu:set_debug(dbgWriter)
	trafoSmoother = AssembledTransformingSmoother2dCPU1(transformedDisc, ilu, rightTrafoDisc)

	return trafoSmoother
end 



--[[
 Schur complement transforming smother
--]]

function CreateBlockTrafoSmoother0(approxSpace, originalDisc, beta)

-- assemble mass matrix (identity)
local schurComplementDisc = DomainDiscretization(approxSpace)
local discP = ConvectionDiffusion("p", "Inner", type)
discP:set_diffusion(0.0)
discP:set_reaction_rate(beta)
schurComplementDisc:add(discP)
--schurComplementDisc:add(BndDisc)
--schurComplementDisc:add(FixPressureDisc)
local trafoSmoother = BlockTransformingSmoother2dCPU1(originalDisc, schurComplementDisc, ILU(), Jacobi(0.66))
trafoSmoother:set_type(0);
return trafoSmoother
end 

function CreateBlockTrafoSmoother1(approxSpace, originalDisc, beta)

-- assemble Laplacian
local schurComplementDisc = DomainDiscretization(approxSpace)
local discP = ConvectionDiffusion("p", "Inner", type)
discP:set_diffusion(beta)
discP:set_reaction_rate(0.0)
schurComplementDisc:add(discP)
--schurComplementDisc:add(BndDisc)
--schurComplementDisc:add(FixPressureDisc)
local trafoSmoother = BlockTransformingSmoother2dCPU1(originalDisc, schurComplementDisc, ILUT())
trafoSmoother:set_type(1);
return trafoSmoother
end 


