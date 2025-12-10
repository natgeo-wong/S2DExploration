### A Pluto.jl notebook ###
# v0.20.19

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
	using Dates, DelimitedFiles
	using StatsBase, MultivariateStats
	using ERA5Reanalysis
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())

	include(srcdir("smoothing.jl"))
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 03a. Vertical Mode Coefficients of ERA5 Pressure Velocity
"

# ╔═╡ 370903b1-1c85-4337-89c0-fd7cc358359d
TableOfContents()

# ╔═╡ 2205ab68-6602-4222-b3d5-1ebf2a7b8082
md"
### A. Loading Station Information
"

# ╔═╡ 8a43991e-cf8e-4304-9d0e-f235a0da1f36
@bind armsite Select([
	"SGP" => "(SGP) Southern Great Plains",
	"BNF" => "(BNF) Bankhead National Forest",
])

# ╔═╡ e282e373-f8a9-4a11-b1a7-e1d270b14090
md"
### B. Defining ERA5 Datasets and Variables
"

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ 73d960be-cbd9-4705-b68d-64bca44478ae
eω = PressureVariable("w")

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Reconstruction of $\omega$ from Vertical Mode Coefficients
"

# ╔═╡ 483383ec-7356-4ad0-a5fb-abb8bf982155
begin
	dω  = read_climatology("SGP",e5ds,eω)
	dt  = dω["valid_time"][:]
	lvl = dω["pressures"][:]
	ωSGP = dω[eω.ncID][:,:]
	close(dω)
	dω  = read_climatology("BNF",e5ds,eω)
	ωBNF = dω[eω.ncID][:,:]
	close(dω)
	md"Loading Surface and Tropopause Pressures"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2001,6,3)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2001,6,8)

# ╔═╡ 52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
md"
### D. What Happens if the Data is Smoothed?
"

# ╔═╡ f1208a2f-47ee-42be-9791-0389a22b9b0d
md"
### E. Decomposition of $\omega$ into Principal Components
"

# ╔═╡ 0ac5abd0-00cf-459b-abac-673717cd1b73
begin
	M = fit(PCA, ωSGP; maxoutdim=2)
	evecSGP = eigvecs(M)
	Y = predict(M,ωSGP)
	rωSGP = reconstruct(M,Y)
	M = fit(PCA, ωBNF; maxoutdim=2)
	evecBNF = eigvecs(M)
	Y = predict(M,ωBNF)
	rωBNF = reconstruct(M,Y)
end

# ╔═╡ 15d830de-0fc7-4d09-8cdd-0081d8b17e0b
begin
	ii = (dt .>= dtbeg) .& (dt .<= dtend)
	f3 = Figure()
	
	ax1 = Axis(
		f3[1,1],width=300,height=200,
		xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	c1 = heatmap!(ax1,1:sum(ii),lvl,ωSGP[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax1,1000,10)
	
	ax2 = Axis(
		f3[1,2],width=300,height=200,
		xticklabelsvisible=false,
		yticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax2,1:sum(ii),lvl,rωSGP[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax2,1000,10)

	ax3 = Axis(
		f3[1,3],width=50,height=200,
		xticklabelsvisible=false,
		yticklabelsvisible=false,xticks=[0]
		# yscale=log10,yticks=[0.1]
	)

	lines!(ax3,evecSGP[:,1]*0.5,lvl)
	lines!(ax3,evecSGP[:,2]*0.25,lvl)
	xlims!(ax3,-0.2,0.2)
	ylims!(ax3,1000,10)

	ax4 = Axis(
		f3[2,1],width=300,height=200,
		xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax4,1:sum(ii),lvl,ωBNF[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax4,1000,10)
	
	ax5 = Axis(
		f3[2,2],width=300,height=200,
		yticklabelsvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax5,1:sum(ii),lvl,rωBNF[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax5,1000,10)

	ax6 = Axis(
		f3[2,3],width=50,height=200,
		yticklabelsvisible=false,xticks=[0]
		# yscale=log10,yticks=[0.1]
	)

	lines!(ax6,evecBNF[:,1]*0.5,lvl)
	lines!(ax6,evecBNF[:,2]*0.25,lvl)
	xlims!(ax6,-0.2,0.2)
	ylims!(ax6,1000,10)

	Label(f3[1,1],halign=:right,valign=:top,padding=(0,7,0,5),L"$\omega_{SGP}$")
	Label(f3[1,2],halign=:right,valign=:top,padding=(0,7,0,5),L"$\omega_{SGP}$ (PC1 + PC2)")
	Label(f3[2,1],halign=:right,valign=:top,padding=(0,7,0,5),L"$\omega_{BNF}$")
	Label(f3[2,2],halign=:right,valign=:top,padding=(0,7,0,5),L"$\omega_{BNF}$ (PC1 + PC2)")

	Label(f3[:,0],"Pressure / hPa",rotation=pi/2)
	Label(f3[3,1:2],"Date")
	Label(f3[3,3],L"$\omega$")
	
	resize_to_layout!(f3)
	f3
end

# ╔═╡ 2dfd4e50-c9d5-40a6-be3b-d9ec2bf9c4ea
md"
### F. Principal Components of Smoothed Data
"

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─73d960be-cbd9-4705-b68d-64bca44478ae
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╠═483383ec-7356-4ad0-a5fb-abb8bf982155
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╟─52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
# ╟─f1208a2f-47ee-42be-9791-0389a22b9b0d
# ╟─0ac5abd0-00cf-459b-abac-673717cd1b73
# ╠═15d830de-0fc7-4d09-8cdd-0081d8b17e0b
# ╟─2dfd4e50-c9d5-40a6-be3b-d9ec2bf9c4ea
