### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5692ccd8-6056-11f0-07da-d1ae989cdac1
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ bd5e50c0-5d29-410e-8389-beb8b636307d
begin
	@quickactivate "S2DExploration"
	using PlutoUI
	using Dates, DelimitedFiles, StatsBase
	using GeoRegions, ERA5Reanalysis
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())

	include(srcdir("smoothing.jl"))
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 03e. Exploration of ERA5 Vertical Shear
"

# ╔═╡ 370903b1-1c85-4337-89c0-fd7cc358359d
TableOfContents()

# ╔═╡ 2205ab68-6602-4222-b3d5-1ebf2a7b8082
md"
### A. Loading Station Information
"

# ╔═╡ 8a43991e-cf8e-4304-9d0e-f235a0da1f36
@bind armsite Select([
	# "SGP" => "(SGP) Southern Great Plains",
	"BNF" => "(BNF) Bankhead National Forest",
])

# ╔═╡ e282e373-f8a9-4a11-b1a7-e1d270b14090
md"
### B. Defining ERA5 Datasets and Variables
"

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
uvar = PressureVariable("u")

# ╔═╡ b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
vvar = PressureVariable("v")

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Timeseries Data for Stations
"

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
	uds = read_climatology("SGP",e5ds,uvar,days=30)
	dtv = uds["valid_time"][:]
	u_SGP = uds[uvar.ncID][:,:]
	close(uds)
	vds = read_climatology("SGP",e5ds,vvar,days=30)
	v_SGP = vds[vvar.ncID][:,:]
	close(vds)
	uds = read_climatology("BNF",e5ds,uvar,days=30)
	u_BNF = uds[uvar.ncID][:,:]
	close(uds)
	vds = read_climatology("BNF",e5ds,vvar,days=30)
	v_BNF = vds[vvar.ncID][:,:]
	close(vds)
	md"Loading 30-day Smoothed Horizontal Velocities (Zonal and Meridional)"
end

# ╔═╡ 5a08ac83-e44f-4f4c-a5c0-bd02ac30a501
md"Year: $(@bind yr PlutoUI.Slider(1980:2024))"

# ╔═╡ ec36729f-5f9a-41b9-84d7-22beaa53f797
md"Month: $(@bind mo PlutoUI.Slider(1:12))"

# ╔═╡ 82ba5e54-c127-42a1-81d4-b3c1bc03d485
md"Day of Month: $(@bind dy PlutoUI.Slider(1:daysinmonth(Date(yr,mo))))"

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dt = DateTime(yr,mo,dy)

# ╔═╡ 6a413bc4-4bb0-48f5-98e0-ffb98430b9e8
md"Do Animation: $(@bind doanim PlutoUI.Slider(0:1))"

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
	# idt = findfirst(dtv .== dt)
	f1 = Figure()

	if isone(doanim)
		for iyr = 1980 : 2024, imo = 1 : 12, idy = 5 : 15 : daysinmonth(Date(iyr,imo))
			idt = findfirst(dtv .== DateTime(iyr,imo,idy))
			dtstr = Dates.format(DateTime(iyr,imo,idy),dateformat"yyyymmdd")
			empty!(f1)
			ax1_1 = Axis(
				f1[1,1],width=400,height=400,limits=(-75,75,-75,75),
				xlabel=L"$u$ / m s$^{-1}$",ylabel=L"$v$ / m s$^{-1}$",
				title = dtstr
			)
		
			for ipre = 12 : 37
				lines!(ax1_1,
					[u_SGP[ipre-1,idt],u_SGP[ipre,idt]],
					[v_SGP[ipre-1,idt],v_SGP[ipre,idt]]
				)
			end
			lines!(ax1_1,[u_SGP[37,idt],0],[v_SGP[37,idt],0])
			
			for ipre = 12 : 37
				lines!(ax1_1,
					[u_BNF[ipre-1,idt],u_BNF[ipre,idt]],
					[v_BNF[ipre-1,idt],v_BNF[ipre,idt]]
				)
			end
			lines!(ax1_1,[u_BNF[37,idt],0],[v_BNF[37,idt],0])
		
			resize_to_layout!(f1)
			CairoMakie.save(plotsdir("03e-exploreverticalshear","$dtstr.png"),f1)
		end
	end
	# f1
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╟─b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─5a08ac83-e44f-4f4c-a5c0-bd02ac30a501
# ╠═ec36729f-5f9a-41b9-84d7-22beaa53f797
# ╟─82ba5e54-c127-42a1-81d4-b3c1bc03d485
# ╟─7cf5b637-8443-45b3-91db-781f01f75810
# ╟─6a413bc4-4bb0-48f5-98e0-ffb98430b9e8
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
