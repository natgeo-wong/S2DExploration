### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 16c7e940-f1b5-11ef-027b-ffed5bebb310
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ 2dd9c606-8ab2-484e-8155-f2c908418fe0
begin
	@quickactivate "S2DExploration"
	using Dates
	using DelimitedFiles
	using ERA5Reanalysis
	using ETOPO
	using GeoRegions
	using PyCall, LaTeXStrings
	using ImageShow, PNGFiles

	pplt = pyimport("proplot")
end

# ╔═╡ 2c95fc06-b3d1-4eef-9d2d-55855a0d3cf6
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	clon  = coast[:,1]
	clat  = coast[:,2]
	ii = (clon .> -70)  .| (clon .< -130) .| (clat .> 47) .| (clat .< 22)
	clon[ii] .= NaN
	clat[ii] .= NaN
	md"Preloading coastline data"
end

# ╔═╡ e7940f44-e357-4f05-9bb2-2b9e0434d1a6
sgp_info = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[:,5:6]

# ╔═╡ 62d834cc-6d05-43b6-9218-7295ac9625b5
bnf_info = readdlm(datadir("ARM_BNF.csv"),',',skipstart=1)[:,5:6]

# ╔═╡ 74a03763-63d0-4a16-9586-f55a5a3aed6f
geo_SGP = GeoRegion("KS_OK",path=srcdir())

# ╔═╡ 8cbbfd33-5093-4720-bd18-bf10372fd357
geo_BNF = GeoRegion("TN_MS_AL",path=srcdir())

# ╔═╡ 7038130a-cfcc-4e3b-a7e7-64157f0ab69b
begin
	lon_SGP,lat_SGP = coordinates(geo_SGP)
	lon_BNF,lat_BNF = coordinates(geo_BNF)
	md"Retrieving GeoRegion coordinates"
end

# ╔═╡ 9d4a687f-43d3-4ef2-b23d-5cd4b1516355
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ 56c7bec1-37f4-4c25-8606-850dc0a367cc
elsd_SGP = getLandSea(e5ds,ERA5Region(geo_SGP))

# ╔═╡ 0ef24a86-627a-4ee3-ab46-793e8ccd3c7f
elsd_BNF = getLandSea(e5ds,ERA5Region(geo_BNF))

# ╔═╡ 8bfec88b-7064-4fb4-908a-d105c71299be
geo = GeoRegion([-105,-80,-80,-105,-105],[28,28,42,42,28])

# ╔═╡ cd72582b-11cf-4697-887c-0b76fd88fe86
etpd = ETOPODataset(path=datadir())

# ╔═╡ 9e7eb0df-58b6-4bad-9e85-709647576360
lsd = getLandSea(etpd,geo)

# ╔═╡ b33690e2-9200-4a00-b010-ca0fe45bfd49
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=25/14,axwidth=3.5)

	axs[1].scatter(sgp_info[:,1],sgp_info[:,2],c="r",s=2,zorder=2)
	axs[1].scatter(bnf_info[:,1],bnf_info[:,2],c="r",s=2,zorder=2)
	
	axs[1].plot(lon_SGP,lat_SGP,c="k",lw=0.5)
	axs[1].plot(lon_BNF,lat_BNF,c="k",lw=0.5)

	axs[1].pcolormesh(
		lsd.lon[1:3:end],lsd.lat[1:3:end],lsd.z[1:3:end,1:3:end]'./1000,
		levels=-1.5:0.1:1.5,extend="both",cmap="bukavu",a=0.6
	)
	c = 
	axs[1].pcolormesh(elsd_SGP.lon,elsd_SGP.lat,elsd_SGP.z'./1000,levels=-1.5:0.1:1.5,extend="both",cmap="bukavu")
	axs[1].pcolormesh(elsd_BNF.lon,elsd_BNF.lat,elsd_BNF.z'./1000,levels=-1.5:0.1:1.5,extend="both",cmap="bukavu")
	
	axs[1].format(
		xlim=(-105,-80),xlabel=L"Longitude / $\degree$",
		ylim=(28,42),ylabel=L"Latitude / $\degree$",ylocator=28:2:42
	)

	fig.colorbar(c,locator=-1.5:0.5:1.5,minorlocator=-1.5:0.1:1.5,label="Topography / km")
	fig.savefig(plotsdir("01-georegions.png"),transparent=false,dpi=300)
	PNGFiles.load(plotsdir("01-georegions.png"))
end

# ╔═╡ Cell order:
# ╟─16c7e940-f1b5-11ef-027b-ffed5bebb310
# ╟─2dd9c606-8ab2-484e-8155-f2c908418fe0
# ╟─2c95fc06-b3d1-4eef-9d2d-55855a0d3cf6
# ╟─e7940f44-e357-4f05-9bb2-2b9e0434d1a6
# ╟─62d834cc-6d05-43b6-9218-7295ac9625b5
# ╟─74a03763-63d0-4a16-9586-f55a5a3aed6f
# ╟─8cbbfd33-5093-4720-bd18-bf10372fd357
# ╟─7038130a-cfcc-4e3b-a7e7-64157f0ab69b
# ╟─9d4a687f-43d3-4ef2-b23d-5cd4b1516355
# ╠═56c7bec1-37f4-4c25-8606-850dc0a367cc
# ╟─0ef24a86-627a-4ee3-ab46-793e8ccd3c7f
# ╟─8bfec88b-7064-4fb4-908a-d105c71299be
# ╟─cd72582b-11cf-4697-887c-0b76fd88fe86
# ╟─9e7eb0df-58b6-4bad-9e85-709647576360
# ╟─b33690e2-9200-4a00-b010-ca0fe45bfd49
