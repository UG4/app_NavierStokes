-- Persistent Data
local multiRefObjects = {
{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};
} -- multiRefObjects
multiRefObjects[15]["y"] = "$\\norm{v_h - v_{{h/2}} }_{ L_2}$";
multiRefObjects[15]["x"] = "h (Gitterweite)";
multiRefObjects[4]["y"] = "$\\norm{v_h - v_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[4]["x"] = "Anzahl Unbekannte";
multiRefObjects[6]["y"] = "$\\norm{p_h - p_{{h/2}} }_{ L_2}$";
multiRefObjects[6]["x"] = "Anzahl Unbekannte";
multiRefObjects[13]["y"] = "$\\norm{p_h - p_{{h/2}} }_{ H^1}$";
multiRefObjects[13]["x"] = "h (Gitterweite)";
multiRefObjects[23]["y"] = "$\\norm{u_h - u_{{h/2}} }_{ H^1}$";
multiRefObjects[23]["x"] = "Anzahl Unbekannte";
multiRefObjects[11]["y"] = "$\\norm{p_h - p_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[11]["x"] = "Anzahl Unbekannte";
multiRefObjects[1]["y"] = "$\\norm{v_h - v_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[1]["x"] = "Anzahl Unbekannte";
multiRefObjects[5]["y"] = "$\\norm{u_h - u_{{h/2}} }_{ L_2}$";
multiRefObjects[5]["x"] = "Anzahl Unbekannte";
multiRefObjects[2]["y"] = "$\\norm{p_h - p_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[2]["x"] = "Anzahl Unbekannte";
multiRefObjects[18]["y"] = "$\\norm{u_h - u_{{h/2}} }_{ L_2}$";
multiRefObjects[18]["x"] = "h (Gitterweite)";
multiRefObjects[12]["y"] = "$\\norm{v_h - v_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[12]["x"] = "h (Gitterweite)";
multiRefObjects[8]["y"] = "$\\norm{p_h - p_{{h/2}} }_{ H^1}$";
multiRefObjects[8]["x"] = "Anzahl Unbekannte";
multiRefObjects[22]["y"] = "$\\norm{p_h - p_{{h/2}} }_{ L_2}$";
multiRefObjects[22]["x"] = "h (Gitterweite)";
multiRefObjects[20]["y"] = "$\\norm{p_h - p_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[20]["x"] = "h (Gitterweite)";
multiRefObjects[19]["y"] = "$\\norm{v_h - v_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[19]["x"] = "h (Gitterweite)";
multiRefObjects[17]["y"] = "$\\norm{u_h - u_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[17]["x"] = "h (Gitterweite)";
multiRefObjects[14]["y"] = "$\\norm{u_h - u_{{h/2}} }_{ H^1}$";
multiRefObjects[14]["x"] = "h (Gitterweite)";
multiRefObjects[21]["y"] = "$\\norm{p_h - p_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[21]["x"] = "h (Gitterweite)";
multiRefObjects[10]["y"] = "$\\norm{v_h - v_{{h/2}} }_{ H^1}$";
multiRefObjects[10]["x"] = "h (Gitterweite)";
multiRefObjects[7]["y"] = "$\\norm{u_h - u_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[7]["x"] = "Anzahl Unbekannte";
multiRefObjects[16]["y"] = "$\\norm{u_h - u_{h_{\\text{min}}} }_{ L_2}$";
multiRefObjects[16]["x"] = "h (Gitterweite)";
multiRefObjects[24]["y"] = "$\\norm{v_h - v_{{h/2}} }_{ L_2}$";
multiRefObjects[24]["x"] = "Anzahl Unbekannte";
multiRefObjects[3]["y"] = "$\\norm{v_h - v_{{h/2}} }_{ H^1}$";
multiRefObjects[3]["x"] = "Anzahl Unbekannte";
multiRefObjects[9]["y"] = "$\\norm{u_h - u_{h_{\\text{min}}} }_{ H^1}$";
multiRefObjects[9]["x"] = "Anzahl Unbekannte";
local obj1 = {
	["plots/dataset/fv1_1_u_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[7];
	};
	["plots/multi/u_fv1_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[16];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[17];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/p_fv1_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[6];
	};
	["plots/multi/p_l-prev_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[22];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/v_fv1_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[4];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[1];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/all_fv1_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[7];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[2];
		};
		[3] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[4];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_all_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[17];
	};
	["plots/dataset/fv1_1_u_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[16];
	};
	["plots/disc/v_fv1_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[4];
	};
	["plots/disc/p_fv1_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[21];
	};
	["plots/disc/v_all_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[19];
	};
	["plots/dataset/fv1_1_p_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[21];
	};
	["plots/disc/u_all_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[7];
	};
	["plots/multi/p_fv1_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[20];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[13];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/v_l-lmax_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[19];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/u_fv1_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[18];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[14];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/dataset/fv1_1_p_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[11];
	};
	["plots/disc/p_all_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[11];
	};
	["plots/multi/u_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[5];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/v_fv1_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[19];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[12];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_fv1_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[16];
	};
	["plots/multi/v_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[3];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/p_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[6];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/p_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[8];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/v_all_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[12];
	};
	["plots/disc/v_all_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[4];
	};
	["plots/dataset/fv1_1_p_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[22];
	};
	["plots/multi/u_fv1_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[7];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[5];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_l-prev_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[14];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/p_all_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[21];
	};
	["plots/disc/v_all_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[10];
	};
	["plots/multi/p_fv1_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[2];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[11];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/dataset/fv1_1_v_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[15];
	};
	["plots/multi/v_fv1_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[12];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[10];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/p_all_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[2];
	};
	["plots/multi/p_l-lmax_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[20];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/u_fv1_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[7];
	};
	["plots/disc/v_fv1_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[10];
	};
	["plots/multi/u_l-lmax_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[17];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/p_fv1_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[21];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[22];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_fv1_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[5];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[23];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/p_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[20];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[21];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/disc/p_fv1_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[13];
	};
	["plots/multi/v_l-lmax_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[12];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/u_l-lmax_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[16];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/u_all_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[16];
	};
	["plots/multi/v_fv1_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[24];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[3];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/v_fv1_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[24];
	};
	["plots/disc/u_fv1_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[17];
	};
	["plots/multi/v_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[4];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/p_fv1_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[8];
	};
	["plots/disc/p_all_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[6];
	};
	["plots/dataset/fv1_1_p_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[13];
	};
	["plots/multi/v_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[1];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[4];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/disc/p_all_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[13];
	};
	["plots/dataset/fv1_1_u_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[14];
	};
	["plots/multi/v_fv1_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[19];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[15];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/all_fv1_l-prev_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[14];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[13];
		};
		[3] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[10];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[9];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[7];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/multi/p_fv1_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[6];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[8];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/v_l-prev_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[15];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/u_all_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[9];
	};
	["plots/multi/p_l-prev_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[13];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/p_fv1_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[22];
	};
	["plots/multi/p_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[11];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/v_fv1_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[4];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[24];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_fv1_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[14];
	};
	["plots/multi/p_fv1_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[2];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[6];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/all_fv1_l-lmax_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[17];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[20];
		};
		[3] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[12];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[17];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[16];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/multi/u_fv1_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[16];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[18];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/p_all_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[20];
	};
	["plots/dataset/fv1_1_v_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[19];
	};
	["plots/multi/all_fv1_l-lmax_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[16];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[21];
		};
		[3] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[19];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/p_all_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[22];
	};
	["plots/disc/v_all_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[1];
	};
	["plots/multi/u_l-prev_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[18];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/dataset/fv1_1_v_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[12];
	};
	["plots/multi/u_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[9];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/dataset/fv1_1_u_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[18];
	};
	["plots/multi/all_fv1_l-prev_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[18];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[22];
		};
		[3] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[15];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_fv1_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[9];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[23];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_fv1_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[9];
	};
	["plots/multi/all_fv1_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[5];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[6];
		};
		[3] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[24];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/v_all_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[15];
	};
	["plots/multi/all_fv1_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[23];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[8];
		};
		[3] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[3];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/dataset/fv1_1_v_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[10];
	};
	["plots/dataset/fv1_1_p_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[20];
	};
	["plots/multi/v_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[10];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[15];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/disc/v_fv1_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[1];
	};
	["plots/disc/u_all_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[18];
	};
	["plots/disc/v_all_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[24];
	};
	["plots/multi/v_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[3];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[24];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/multi/u_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[23];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/p_fv1_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[20];
	};
	["plots/multi/p_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[13];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[22];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/dataset/fv1_1_v_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[3];
	};
	["plots/multi/u_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[23];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[5];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/disc/p_fv1_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[11];
	};
	["plots/multi/u_fv1_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[7];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[9];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/p_l-lmax_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[11];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[2];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/disc/v_fv1_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[12];
	};
	["plots/dataset/fv1_1_p_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[2];
	};
	["plots/multi/v_l-prev_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[10];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/p_all_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[8];
	};
	["plots/disc/u_fv1_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[18];
	};
	["plots/multi/p_fv1_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[22];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[13];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_all_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[23];
	};
	["plots/multi/v_fv1_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[15];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[10];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/v_all_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[3];
	};
	["plots/multi/p_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[2];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/v_fv1_l-prev_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[15];
	};
	["plots/multi/v_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[24];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/disc/u_all_l-prev_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[14];
	};
	["plots/dataset/fv1_1_v_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[24];
	};
	["plots/disc/v_fv1_l-lmax_l2_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[19];
	};
	["plots/dataset/fv1_1_p_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[6];
	};
	["plots/dataset/fv1_1_p_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[8];
	};
	["plots/disc/p_fv1_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[2];
	};
	["plots/multi/p_l-lmax_l2_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[21];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/dataset/fv1_1_v_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[4];
	};
	["plots/multi/v_fv1_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[1];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[3];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/dataset/fv1_1_u_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[9];
	};
	["plots/dataset/fv1_1_u_l-lmax_h1_h"] = {
		[1] = {
			[1] = 2;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[17];
	};
	["plots/disc/u_all_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[5];
	};
	["plots/multi/u_l-prev_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[14];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[18];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/multi/v_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[12];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[19];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/multi/p_fv1_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[11];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[8];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_fv1_h1_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[17];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[14];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/u_l-lmax_l2_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[7];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/multi/p_fv1_l-lmax_h"] = {
		[1] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[21];
		};
		[2] = {
			[1] = {
				[1] = 2;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[20];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/multi/all_fv1_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_u_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[9];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[11];
		};
		[3] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[1];
		};
		["multiplot"] = {
			["cols"] = 1;
		};
	};
	["plots/disc/u_fv1_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[23];
	};
	["plots/disc/v_fv1_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[3];
	};
	["plots/multi/p_l-prev_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[8];
		};
		[2] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_p_l-prev_l2.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[6];
		};
		["multiplot"] = {
			["cols"] = 2;
		};
	};
	["plots/dataset/fv1_1_u_l-prev_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[23];
	};
	["plots/dataset/fv1_1_v_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[1];
	};
	["plots/multi/v_l-lmax_h1_DoFs"] = {
		[1] = {
			[1] = {
				[1] = 1;
				[2] = 3;
				["file"] = "data/error_fv1_1_v_l-lmax_h1.dat";
				["style"] = "linespoints";
				["label"] = "fv1 $\\mathbb{P}_{1}$";
			};
			["label"] = multiRefObjects[1];
		};
		["multiplot"] = {
			["rows"] = 1;
		};
	};
	["plots/dataset/fv1_1_u_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[5];
	};
	["plots/disc/u_fv1_l-prev_l2_DoFs"] = {
		[1] = {
			[1] = 1;
			[2] = 3;
			["file"] = "data/error_fv1_1_u_l-prev_l2.dat";
			["style"] = "linespoints";
			["label"] = "fv1 $\\mathbb{P}_{1}$";
		};
		["label"] = multiRefObjects[5];
	};
}
return obj1
