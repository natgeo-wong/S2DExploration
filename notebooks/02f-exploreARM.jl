### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 5692ccd8-6056-11f0-07da-d1ae989cdac1
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ bd5e50c0-5d29-410e-8389-beb8b636307d
begin
	@quickactivate "S2DExploration"
	using Dates
	using DelimitedFiles
	using NASAMergedTb
	using NCDatasets
	using Statistics
	
	using CairoMakie
	using LaTeXStrings
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02. Explore Some Data
"

# ╔═╡ 1a2cd78c-4cdb-4ee6-8f69-f649fbdb535c
begin
	slon,slat = readdlm(datadir("ARM_SGP.csv"),',',skipstart=1)[1,5:6]
	md"Loading station coordinates ..."
end

# ╔═╡ 445d35d6-1846-4e96-ba6d-0448b25abbdc
btd = TbDataset(start=Date(2001),stop=Date(2024),path=datadir())

# ╔═╡ 4d7f159b-9df3-4b41-8c42-eef810f4ddea
geo = GeoRegion("KS_OK",path=srcdir())

# ╔═╡ 1b9fedd0-4127-461b-80a0-a8434ddd26ca
begin
	ds  = read(btd,geo,Date(2011,6,1))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	tb  = ds["Tb"][:,:,:]
	close(ds)
end

# ╔═╡ bb65b62e-33d0-4f33-b12a-dbf5d7178adc
ilon = argmin(abs.(lon.-slon))

# ╔═╡ 7bc627ca-b580-4ded-8716-18ff156b7b56
ilat = argmin(abs.(lat.-slat))

# ╔═╡ 215c442a-634a-472e-b0bb-61ad561735bb
itb = tb[ilon,ilat,:]

# ╔═╡ e0e3cea5-c898-43d3-99b9-8860efc3bd48
tvec = DateTime(2011,6,1,0,0,0) : Minute(30) : DateTime(2011,6,1,23,59,0)

# ╔═╡ f2715353-7ff5-44e0-9eaa-e929aa10b4d1
begin
	ads = NCDataset("/scratch/an4462/videos/2011/sgpdlfptC1.b1.20110602.050001.cdf")
	t   = ads["time"][:]
	z   = ads["range"][:]
	w   = nomissing(ads["radial_velocity"][:,:],NaN)
	close(ads)
end

# ╔═╡ c3f9c5eb-b8a3-4bef-b69e-a4333b5d9fc0
begin
	fig = Figure()
	
	ax1 = Axis(
	    fig[1,1],width=500,height=250,
	    # limits=(geo.W-0.5,geo.E+0.5,geo.S-0.5,geo.N+0.5)
	)
	c = heatmap!(ax1,w',colorrange=(-2,2),colormap=:RdBu)

	# Colorbar(fig[1,2], c)
	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╠═1a2cd78c-4cdb-4ee6-8f69-f649fbdb535c
# ╠═445d35d6-1846-4e96-ba6d-0448b25abbdc
# ╠═4d7f159b-9df3-4b41-8c42-eef810f4ddea
# ╠═1b9fedd0-4127-461b-80a0-a8434ddd26ca
# ╠═bb65b62e-33d0-4f33-b12a-dbf5d7178adc
# ╠═7bc627ca-b580-4ded-8716-18ff156b7b56
# ╠═215c442a-634a-472e-b0bb-61ad561735bb
# ╠═e0e3cea5-c898-43d3-99b9-8860efc3bd48
# ╠═f2715353-7ff5-44e0-9eaa-e929aa10b4d1
# ╠═c3f9c5eb-b8a3-4bef-b69e-a4333b5d9fc0
