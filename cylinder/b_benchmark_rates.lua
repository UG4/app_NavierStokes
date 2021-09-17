if bBenchmarkRates then

	local dom = CreateDomain()
	local plots = {}	

	local function ComputeSpace(discType, p, ppress, minLev, maxLev)

		local file = table.concat({"dc",discType,p,ppress},"_")..".dat"
		--[[
		local discLabel = discType.." Q_"..p.."/Q_"..(ppress)
		local CDLabel = "|C_D - C_D_h|"
		local CLLabel = "|C_L - C_L_h|"
		local DeltaPLabel = "|P - P_L_h|"
		--]]
		local discLabel = discType.." $\\mathbb{Q}_{"..p.."}/\\mathbb{Q}_{"..(ppress).."}$"
		local CDLabel = "$|c_{w,h} - c_w^{\\text{Ref}}|$"
		local CLLabel = "$|c_{a,h} - c_a^{\\text{Ref}}|$"
		local DeltaPLabel = "$|\\Delta_{p,h} - \\Delta_p^{\\text{Ref}}|$"
		
		local function addPlot(name, dataset, label)
			plots[name] = plots[name] or {}
			table.insert( plots[name], dataset)			
			plots[name].label = label			
		end
		
		addPlot("CD_DoF", {label=discLabel, file=file, style="linespoints", 1, 3},
				{ x = "Anzahl Unbekannte", y = CDLabel})

		addPlot("CL_DoF", {label=discLabel, file=file, style="linespoints", 1, 4},
				{ x = "Anzahl Unbekannte", y = CLLabel})

		addPlot("DeltaP_DoF", {label=discLabel, file=file, style="linespoints", 1, 5},
				{ x = "Anzahl Unbekannte", y = DeltaPLabel})

		addPlot("CD_h", {label=discLabel, file=file, style="linespoints", 2, 3},
				{ x = "h (Gitterweite)", y = CDLabel})

		addPlot("CL_h", {label=discLabel, file=file, style="linespoints", 2, 4},
				{ x = "h (Gitterweite)", y = CLLabel})

		addPlot("DeltaP_h", {label=discLabel, file=file, style="linespoints", 2, 5},
				{ x = "h (Gitterweite)", y = DeltaPLabel})

		if not util.HasParamOption("-replot") then

			local approxSpace = util.ns.CreateApproxSpace(dom, discType, p, ppress)
			local domainDisc = CreateDomainDisc(approxSpace, discType, p, ppress)
			local solver = CreateSolver(approxSpace, discType, p)
		
			local h, DoFs, level = {}, {}, {}	
			local meas = {CD = {}, CL = {}, DeltaP = {}}
			for n, _ in pairs(meas) do
				meas[n] = {value = {}, prev = {error = {}, rate = {}}, exact = {error = {}, rate = {}}}
			end	
			
			local uPrev = nil
			for lev = minLev, maxLev do
				write("\n>> Computing Level "..lev..", "..discType..", "..p..", "..ppress..".\n")
		
				local u = GridFunction(approxSpace, lev)

				if uPrev ~= nil and discType ~= "fvcr" and discType ~= "fecr" then	
					Prolongate(u, uPrev);
					AdjustMeanValue(u, "p")
					u:enforce_consistent_type()
					u:check_storage_type()
				else
					u:set(0)	
				end
				
				if lev > minLev + 2 then globalNSDisc:set_exact_jacobian(true) end
				
				ComputeNonLinearSolution(u, domainDisc, solver)
				u:check_storage_type()

				local FctCmp = approxSpace:names()
				local VelCmp = {}
				
				for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
				vtkWriter = VTKOutput()
				vtkWriter:select(VelCmp, "velocity")
				vtkWriter:select("u", "u")
				vtkWriter:select("v", "v")
				vtkWriter:select("p", "pressure")
				vtkWriter:print("Cylinder"..table.concat({discType,p,ppress,"lev",lev},"_"), u)
		
				-- h/DoF statistic
				DoFs[lev] = u:num_dofs()
				h[lev] =  MaxElementDiameter(dom, lev) 
				level[lev] = lev		
				
				-- C_D / C_L
				local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p*p+5)
				meas.CD.value[lev] = 2*DL[1]/(Umean2*L)
				meas.CL.value[lev] = 2*DL[2]/(Umean2*L)
			
				-- Delta_P
				local PEval = GlobalGridFunctionNumberData(u, "p")
				meas.DeltaP.value[lev] = PEval:evaluate_global({0.15, 0.2}) - PEval:evaluate_global( {0.25, 0.2} )
		
				for n, _ in pairs(meas) do
					local quant = meas[n]
					quant.exact.error[lev] = math.abs(quant.value[lev] - ref[n])
					
					if lev > minLev then
						quant.prev.error[lev-1] = math.abs(quant.value[lev] - quant.value[lev-1])
					end			
					for _, t in ipairs({"exact", "prev"}) do
						local type = quant[t]
						if type.error[lev-2] ~= nil and type.error[lev-1] ~= nil then
							local fac = type.error[lev-2] / type.error[lev-1]
							type.rate[lev-1] = math.log(fac) / math.log(2) --math.log(h[lev-1]/h[lev])
						end			
						if type.error[lev-1] ~= nil and type.error[lev] ~= nil then
							local fac = type.error[lev-1] / type.error[lev]
							type.rate[lev] = math.log(fac) / math.log(2) --math.log(h[lev-1]/h[lev])
						end			
					end
				end
						
				-- print
				local values = {level, h, DoFs}
				local heading = {"L", "h", "#DoFs"}
				local format = {"%d", "%.2e", "%d"}
		
				for n, _ in pairs(meas) do
					local quant = meas[n]
					table.append(values, {quant.value}) 
					table.append(heading,{n})
					table.append(format, {"%.8f"})
		
					for _, t in ipairs({"exact", "prev"}) do
						local type = quant[t]
						table.append(values, {type.error, type.rate}) 
						table.append(heading,{n.." "..t, "rate"})
						table.append(format, {"%.3e", "%.3f"})
					end
				end
		
				write("\n>> Stats: lev "..lev..", "..discType..", "..p..", "..ppress..".\n")
				table.print(values, {heading = heading, format = format, 
									 hline = true, vline = true, forNil = "--"})
									 
				uPrev = u
				--collectgarbage()
			end
			
			local cols = {DoFs, h, meas.CD.exact.error, meas.CL.exact.error, meas.DeltaP.exact.error}
			gnuplot.write_data(file, cols)	
			
			approxSpace, domainDisc, solver = nil, nil, nil
			collectgarbage()
		end
	end
	
--	ComputeSpace("fecr", 1, 0, numPreRefs, numRefs)	
--	ComputeSpace("fvcr", 1, 0, numPreRefs, numRefs)	
--	ComputeSpace("fe", 1, 1, numPreRefs, numRefs)	
--	ComputeSpace("fe", 2, 1, numPreRefs, numRefs-1)	
--	ComputeSpace("fe", 3, 2, numPreRefs, numRefs-2)	
	ComputeSpace("fv1", 1, 1, numPreRefs, numRefs)	
--	ComputeSpace("fv", 2, 1, numPreRefs, numRefs-1)	
--	ComputeSpace("fv", 3, 2, numPreRefs, numRefs-2)	
	
	if util.HasParamOption("-replot") then
		local texOptions = {	
		
			size = 				{12.5, 6.75}, -- the size of canvas (i.e. plot)
			sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
			font = 				"Arial",
			fontsize =			12,
			fontscale = 		0.7,
			
			logscale = 			true,
			grid = 				"lc rgb 'grey70' lt 0 lw 1", 
			linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
								linewidth = 3, pointsize = 1.3},
			border = 			" back lc rgb 'grey40' lw 2",
			decimalsign = 		",",
			key =	 			"on box lc rgb 'grey40' left bottom Left reverse spacing 2 width 1.1 samplen 2 height 0.5",
			tics =	 			{x = "nomirror out scale 0.75 format '%.te%01T' font ',8'",
								 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
			mtics =				5,
			slope = 			{dy = 3, quantum = 0.5, at = "last"},
			padrange = 			{ x = {0.6, 2}, y = {0.6, 2}},
		}

		local pdfOptions = {	
		
			size = 				{12.5, 9.75}, -- the size of canvas (i.e. plot)
			sizeunit =			"cm", -- one of: cm, mm, {in | inch}, {pt | pixel}
			font = 				"Arial",
			fontsize =			8,
			
			logscale = 			true,
			grid = 				"lc rgb 'grey70' lt 0 lw 1", 
			linestyle =			{colors = gnuplot.RGBbyLinearHUEandLightness(8, 1, 360+40, 85, 0.4, 0.4), 
								linewidth = 3, pointsize = 1.3},
			border = 			" back lc rgb 'grey40' lw 2",
			decimalsign = 		",",
			key =	 			"on box lc rgb 'grey40' left bottom Left reverse spacing 2 width 1.1 samplen 2 height 0.5",
			tics =	 			{x = "nomirror out scale 0.75 format '%g' font ',8'",
								 y = "10 nomirror out scale 0.75 format '%.te%01T' font ',8'"}, 
			mtics =				5,
			slope = 			{dy = 3, quantum = 0.25, at = "last"},
			padrange = 			{ x = {0.6, 10}, y = {0.1, 1.1}},
		}
	
		for name,data in pairs(plots) do
			gnuplot.plot(name..".pdf", data, pdfOptions)
			gnuplot.plot(name..".tex", data, texOptions)
		end
	end
end
