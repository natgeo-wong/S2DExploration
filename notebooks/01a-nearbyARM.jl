### A Pluto.jl notebook ###
# v0.20.17

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

# ╔═╡ 16c7e940-f1b5-11ef-027b-ffed5bebb310
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ 2dd9c606-8ab2-484e-8155-f2c908418fe0
begin
	@quickactivate "S2DExploration"
	using PlutoUI
	using Dates, DelimitedFiles
	using GeoRegions, ERA5Reanalysis, ETOPO
	using Distances
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 2f6fac9c-a075-403b-a696-1fcb095f72d8
md"
# 01a. Neighbouring ARM Stations
"

# ╔═╡ 74a03763-63d0-4a16-9586-f55a5a3aed6f
geo_SGP = GeoRegion("plot_SGP",path=srcdir())

# ╔═╡ 8cbbfd33-5093-4720-bd18-bf10372fd357
geo_BNF = GeoRegion("plot_BNF",path=srcdir())

# ╔═╡ 7038130a-cfcc-4e3b-a7e7-64157f0ab69b
begin
	sgp_beg  = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[:,3]
	sgp_end  = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[:,4]
	sgp_end[sgp_end.=="N/A"] .= "$(Date(now()))"
	sgp_beg  = Date.(sgp_beg); sgp_end  = Date.(sgp_end)
	ii = (sgp_beg .< Date(2021)) .& (sgp_end .> Date(2025))
	
	sgp_info = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[ii,[1,2,5,6]]
	bnf_info = readdlm(datadir("ARM_BNF.csv"),',',skipstart=1)[:,[1,2,5,6]]
	sgp_info = sgp_info[.!isnan.(sgp_info[:,3]),:]
	bnf_info = bnf_info[.!isnan.(bnf_info[:,3]),:]

	nsgp = size(sgp_info,1)
	nbnf = size(bnf_info,1)
	md"Loading ARM station coordinates in the Southern Great Plains"
end

# ╔═╡ 72ac56df-7c31-408a-85c4-7d7a4930676a
begin
	elon,elat = nativelonlat()
	md"Loading the native longitude/latitude points for the ERA5 Reanalysis model ..."
end

# ╔═╡ cd72582b-11cf-4697-887c-0b76fd88fe86
etpd = ETOPODataset(path=datadir())

# ╔═╡ 9e7eb0df-58b6-4bad-9e85-709647576360
lsd_SGP = getLandSea(etpd,geo_SGP)

# ╔═╡ dd2d3404-7cc5-40cf-a961-156598a6c38f
lsd_BNF = getLandSea(etpd,geo_BNF)

# ╔═╡ b08f15d1-9820-4872-97cb-1e186862c08c
@bind iisgp Select([ii => "($(sgp_info[ii,1])) $(sgp_info[ii,2])" for ii in 1 : nsgp])

# ╔═╡ 7a92a478-d3a4-4163-8af9-7272afb88469
@bind iibnf Select([ii => "($(bnf_info[ii,1])) $(bnf_info[ii,2])" for ii in 1 : nbnf])

# ╔═╡ 9ec02d97-fd8b-4af7-b514-5d581786e112
begin
	iilonsgp = sgp_info[iisgp,3]
	iilatsgp = sgp_info[iisgp,4]
	sgp_dist = zeros(nsgp)
	for istn in 1 : nsgp
		sgp_dist[istn] = haversine((iilonsgp,iilatsgp),(sgp_info[istn,3],sgp_info[istn,4])) /1e3
	end
end

# ╔═╡ 89050b8b-3e1f-448b-8aa2-76518e6011df
begin
	iilonbnf = bnf_info[iibnf,3]
	iilatbnf = bnf_info[iibnf,4]
	bnf_dist = zeros(nbnf)
	for istn in 1 : nbnf
		bnf_dist[istn] = haversine((iilonbnf,iilatbnf),(bnf_info[istn,3],bnf_info[istn,4])) /1e3
	end
end

# ╔═╡ b33690e2-9200-4a00-b010-ca0fe45bfd49
begin
	fig = Figure()
	
	ax1 = Axis(
		fig[1,1],width=200,height=200,title="Southern Great Plains",
		limits=(-98.2,-95.8,35.3,37.7),xticks=-100:0.5:95,yticks=34:0.5:39
	)
	heatmap!(ax1,lsd_SGP.lon,lsd_SGP.lat,lsd_SGP.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
	scatter!(ax1,sgp_info[:,3],sgp_info[:,4],color=:black)
	scatter!(ax1,sgp_info[sgp_dist.<20,3],sgp_info[sgp_dist.<20,4],color=:blue)
	scatter!(ax1,sgp_info[iisgp,3],sgp_info[iisgp,4],color=:red,markersize=12)
	scatter!(ax1,elon.-360,elat,markersize=5,color=:grey)
	
	ax2 = Axis(
		fig[1,2],width=200,height=200,title="Bankhead National Forest",
		limits=(-88.2,-85.8,33.3,35.7),xticks=-90:0.5:-85,yticks=32:0.5:37,
		yaxisposition=:right
	)
	c = heatmap!(ax2,lsd_BNF.lon,lsd_BNF.lat,lsd_BNF.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
	scatter!(ax2,bnf_info[:,3],bnf_info[:,4],color=:black)
	scatter!(ax2,bnf_info[bnf_dist.<20,3],bnf_info[bnf_dist.<20,4],color=:blue)
	scatter!(ax2,bnf_info[iibnf,3],bnf_info[iibnf,4],color=:red,markersize=12)
	scatter!(ax2,elon.-360,elat,markersize=5,color=:grey)
	
	Label(fig[2,:],L"Longitude / $\degree$")
	Label(fig[1,0],L"Latitude / $\degree$",rotation=pi/2)

	Colorbar(fig[1,3], c, label="Topographic Height / km")
	resize_to_layout!(fig)
	fig
end

# ╔═╡ 6d8aef7c-e9c2-4b0f-b27a-707cf910b721
CairoMakie.save(plotsdir("01-nearbyARM-$(sgp_info[iisgp,1])$(bnf_info[iibnf,1]).png"),fig)

# ╔═╡ Cell order:
# ╟─2f6fac9c-a075-403b-a696-1fcb095f72d8
# ╟─16c7e940-f1b5-11ef-027b-ffed5bebb310
# ╟─2dd9c606-8ab2-484e-8155-f2c908418fe0
# ╟─74a03763-63d0-4a16-9586-f55a5a3aed6f
# ╟─8cbbfd33-5093-4720-bd18-bf10372fd357
# ╟─7038130a-cfcc-4e3b-a7e7-64157f0ab69b
# ╟─72ac56df-7c31-408a-85c4-7d7a4930676a
# ╟─cd72582b-11cf-4697-887c-0b76fd88fe86
# ╟─9e7eb0df-58b6-4bad-9e85-709647576360
# ╟─dd2d3404-7cc5-40cf-a961-156598a6c38f
# ╟─b08f15d1-9820-4872-97cb-1e186862c08c
# ╟─7a92a478-d3a4-4163-8af9-7272afb88469
# ╟─9ec02d97-fd8b-4af7-b514-5d581786e112
# ╟─89050b8b-3e1f-448b-8aa2-76518e6011df
# ╟─b33690e2-9200-4a00-b010-ca0fe45bfd49
# ╟─6d8aef7c-e9c2-4b0f-b27a-707cf910b721
