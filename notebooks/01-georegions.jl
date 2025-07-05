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
	md"Preloading coastline data"
end

# ╔═╡ e7940f44-e357-4f05-9bb2-2b9e0434d1a6
sgp_info = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[:,5:6]

# ╔═╡ 62d834cc-6d05-43b6-9218-7295ac9625b5
bnf_info = readdlm(datadir("ARM_BNF.csv"),',',skipstart=1)[:,5:6]

# ╔═╡ 74a03763-63d0-4a16-9586-f55a5a3aed6f
if isID("SGP",path=srcdir())
	geo_SGP = GeoRegion("SGP",path=srcdir())
else
	geo_SGP = GeoRegion(
		[-100,-95,-95,-100,-100], [34,34,39,39,34],
		ID = "SGP", pID = "GLB", name = "Southern Great Plains",
		save = true, path = srcdir()
	)
end

# ╔═╡ 41e93beb-458e-4606-aba8-ff1a25917fb4
if isID("BNF",path=srcdir())
	geo_BNF = GeoRegion("BNF",path=srcdir())
else
	geo_BNF = GeoRegion(
		[-87.75,-86.75,-86.75,-87.75,-87.75], [34,34,35,35,34],
		ID = "BNF", pID = "GLB", name = "Bankhead National Forest",
		save = true, path = srcdir()
	)
end

# ╔═╡ 4c2ab4fe-908f-44d1-afaf-0b6aa853bfa6
geo_OK = GeoRegion("OK",path=srcdir())

# ╔═╡ cc516500-cc57-4a35-8ca2-62c1f6d683fd
geo_KS = GeoRegion("KS",path=srcdir())

# ╔═╡ 30bda449-41b6-4dae-afad-40c56ee6cc72
geo_AL = GeoRegion("AL",path=srcdir())

# ╔═╡ 7038130a-cfcc-4e3b-a7e7-64157f0ab69b
begin
	lon_SGP,lat_SGP = coordinates(geo_SGP)
	lon_BNF,lat_BNF = coordinates(geo_BNF)
	lon_OK ,lat_OK  = coordinates(geo_OK)
	lon_KS ,lat_KS  = coordinates(geo_KS)
	lon_AL ,lat_AL  = coordinates(geo_AL)
	md"Retrieving GeoRegion coordinates"
end

# ╔═╡ b33690e2-9200-4a00-b010-ca0fe45bfd49
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=25/12)
	
	axs[1].plot(clon,clat,c="k",lw=0.5)

	axs[1].scatter(sgp_info[:,1],sgp_info[:,2],s=2)
	axs[1].scatter(bnf_info[:,1],bnf_info[:,2],s=2)
	
	axs[1].plot(lon_OK,lat_OK,c="k",lw=1)
	axs[1].plot(lon_KS,lat_KS,c="k",lw=1)
	axs[1].plot(lon_AL,lat_AL,c="k",lw=1)

	axs[1].plot(lon_SGP,lat_SGP)
	axs[1].plot(lon_BNF,lat_BNF)
	
	axs[1].format(
		xlim=(-105,-80),xlabel=L"Longitude / $\degree$",
		ylim=(30,42),ylabel=L"Latitude / $\degree$",ylocator=30:2:42
	)
	
	fig.savefig("test.png",transparent=false,dpi=200)
	PNGFiles.load("test.png")
end

# ╔═╡ Cell order:
# ╟─16c7e940-f1b5-11ef-027b-ffed5bebb310
# ╟─2dd9c606-8ab2-484e-8155-f2c908418fe0
# ╟─2c95fc06-b3d1-4eef-9d2d-55855a0d3cf6
# ╟─e7940f44-e357-4f05-9bb2-2b9e0434d1a6
# ╟─62d834cc-6d05-43b6-9218-7295ac9625b5
# ╟─74a03763-63d0-4a16-9586-f55a5a3aed6f
# ╟─41e93beb-458e-4606-aba8-ff1a25917fb4
# ╟─4c2ab4fe-908f-44d1-afaf-0b6aa853bfa6
# ╟─cc516500-cc57-4a35-8ca2-62c1f6d683fd
# ╟─30bda449-41b6-4dae-afad-40c56ee6cc72
# ╠═7038130a-cfcc-4e3b-a7e7-64157f0ab69b
# ╠═b33690e2-9200-4a00-b010-ca0fe45bfd49
