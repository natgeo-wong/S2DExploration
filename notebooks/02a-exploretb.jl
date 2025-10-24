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
	using GeoRegions, NASAMergedTb
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02a. Exploration of NASA Brightness Temperature Data
"

# ╔═╡ 370903b1-1c85-4337-89c0-fd7cc358359d
TableOfContents()

# ╔═╡ 2205ab68-6602-4222-b3d5-1ebf2a7b8082
md"
### A. Loading Station Information
"

# ╔═╡ baaf05e1-f7a8-49f1-850d-d42cb2ed0f2e
begin	
	sgp_info = readdlm(datadir("ARMstations_SGP.csv"),',',skipstart=1)[:,1:6]
	bnf_info = readdlm(datadir("ARMstations_BNF.csv"),',',skipstart=1)[:,1:6]
	sgp_info = sgp_info[.!isnan.(sgp_info[:,5]),:]
	bnf_info = bnf_info[.!isnan.(bnf_info[:,5]),:]
	
	sgp_info[sgp_info[:,4].=="N/A",4] .= "$(Date(now()))"
	sgp_info[:,3]  = Date.(sgp_info[:,3])
	sgp_info[:,4]  = Date.(sgp_info[:,4])

	nsgp = size(sgp_info,1)
	nbnf = size(bnf_info,1)
	md"Loading ARM station coordinates in the Southern Great Plains"
end

# ╔═╡ 8a43991e-cf8e-4304-9d0e-f235a0da1f36
@bind armsite Select([
	"SGP" => "(SGP) Southern Great Plains",
	"BNF" => "(BNF) Bankhead National Forest",
])

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### B. Preliminary Spatial Exploration
"

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
btd = TbDataset(start=Date(2001),stop=Date(2001),path=datadir())

# ╔═╡ fdf12d12-ce1e-46bf-8ad2-550a615364e3
geo = GeoRegion(armsite,path=srcdir())

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
	ds  = read(btd,geo,Date(2021,6,2))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	tb  = ds["Tb"][:,:,8]
	close(ds)
end

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
	fig = Figure()
	
	ax1 = Axis(
	    fig[1,1],width=250*(geo.E-geo.W+1)/(geo.N-geo.S+1),height=250,
	    limits=(geo.W-0.5,geo.E+0.5,geo.S-0.5,geo.N+0.5)
	)
	c = heatmap!(ax1,lon,lat,tb,colorrange=(200,270),colormap=:Greys)
	if armsite == "SGP"
		scatter!(ax1,sgp_info[:,5],sgp_info[:,6])
		scatter!(ax1,sgp_info[1,5],sgp_info[1,6])
	else
		scatter!(ax1,bnf_info[:,5],bnf_info[:,6])
		scatter!(ax1,bnf_info[1,5],bnf_info[1,6])
	end
	
	Label(fig[2,1],L"Longitude / $\degree$")
	Label(fig[1,0],L"Latitude / $\degree$",rotation=pi/2)

	Colorbar(fig[1,2], c)
	resize_to_layout!(fig)
	fig
end

# ╔═╡ bba26d40-6808-485e-9834-82ea4dc6c176
md"
### C. Timeseries Data for Stations
"

# ╔═╡ 1a48c947-198a-4e85-af2d-87dc2eefe94a
if armsite == "SGP"
	@bind iistn1 Select([ii => "($(sgp_info[ii,1])) $(sgp_info[ii,2])" for ii in 1 : nsgp])
else
	@bind iistn1 Select([ii => "($(bnf_info[ii,1])) $(bnf_info[ii,2])" for ii in 1 : nbnf])
end

# ╔═╡ 44c249c0-0805-438b-a7d7-a5a67a62be18
if armsite == "SGP"
	@bind iistn2 Select([ii => "($(sgp_info[ii,1])) $(sgp_info[ii,2])" for ii in 1 : nsgp])
else
	@bind iistn2 Select([ii => "($(bnf_info[ii,1])) $(bnf_info[ii,2])" for ii in 1 : nbnf])
end

# ╔═╡ 17897c14-ebce-4947-8c85-668d4818423f
if armsite == "SGP"
	stnID1 = sgp_info[iistn1,1];
	stnID2 = sgp_info[iistn2,1];
else
	stnID1 = bnf_info[iistn1,1];
	stnID2 = bnf_info[iistn2,1];
end;

# ╔═╡ d9818f36-6dd5-45d1-b56b-42392ed79049
begin
	tds = NCDataset(datadir(
		"mergedIR","timeseries","$(armsite)_$(stnID1)-Tb-20010101-20241231.nc"
	))
	dt1 = tds["valid_time"][:]
	Tb1 = tds["Tb"][:]
	close(tds)
	tds = NCDataset(datadir(
		"mergedIR","timeseries","$(armsite)_$(stnID2)-Tb-20010101-20241231.nc"
	))
	dt2 = tds["valid_time"][:]
	Tb2 = tds["Tb"][:]
	close(tds)
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2001)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2004)

# ╔═╡ 634c8c07-1aac-4763-86ae-b2b13519164a
begin
	ii = (dt1 .>= dtbeg) .& (dt2 .<= dtend)
	f2 = Figure()

	if armsite == "SGP"
		stn1str = sgp_info[iistn1,1]
		stn2str = sgp_info[iistn2,1]
	else
		stn1str = bnf_info[iistn1,1]
		stn2str = bnf_info[iistn2,1]
	end
	
	ax2_1 = Axis(
		f2[1,1],width=400,height=200,yticks=200:20:320,
		xlabel="Date",ylabel="Brightness Temperature (Tb) / K"
	)

	lines!(ax2_1,dt1[ii],Tb1[ii],label=stn1str,color=:black,linewidth=3)
	lines!(ax2_1,dt2[ii],Tb2[ii],label=stn2str)
	xlims!(ax2_1,dtbeg,dtend)
	ylims!(ax2_1,190,330)

	axislegend()
	
	ax2_2 = Axis(
		f2[1,2],width=200,height=200,limits=(190,330,190,330),
		xticks=200:20:320,yticks=200:20:320,yaxisposition=:right,
		xlabel="$stn1str Tb / K",ylabel="$stn2str Tb / K",
	)

	bins = 190:2:330
	h = fit(Histogram,(Tb1[ii],Tb2[ii]),(bins,bins))
	fmat = Float64.(h.weights)
	fmat[iszero.(fmat)] .= NaN

	c2 = heatmap!(ax2_2,bins,bins,log10.(fmat),colorrange=(0,3))
	lines!(ax2_2,[190,330],[190,330],color=:grey,linewidth=0.5)
	translate!(c2, 0, 0, -100)

	Colorbar(f2[1,3], c2)
	resize_to_layout!(f2)
	f2
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─baaf05e1-f7a8-49f1-850d-d42cb2ed0f2e
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─fdf12d12-ce1e-46bf-8ad2-550a615364e3
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─bba26d40-6808-485e-9834-82ea4dc6c176
# ╟─1a48c947-198a-4e85-af2d-87dc2eefe94a
# ╟─44c249c0-0805-438b-a7d7-a5a67a62be18
# ╟─17897c14-ebce-4947-8c85-668d4818423f
# ╟─d9818f36-6dd5-45d1-b56b-42392ed79049
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╟─634c8c07-1aac-4763-86ae-b2b13519164a
