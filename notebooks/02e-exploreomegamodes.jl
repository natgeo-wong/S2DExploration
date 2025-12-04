### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
eω1 = SingleVariable("ω1",path=srcdir())

# ╔═╡ b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
eω2 = SingleVariable("ω2",path=srcdir())

# ╔═╡ 2a9efd0b-1e82-4c08-8239-451ab6ebcea0
eptrop = SingleVariable("ptrop",path=srcdir())

# ╔═╡ 43502ac5-f6be-4858-8969-4e3b2e9498cd
esp = SingleVariable("sp")

# ╔═╡ 42aa31cf-796b-4b92-8eee-1e53c9d96765
function ωreconstruction(c1,c2,ptrop,sp)

	pre = ERA5Reanalysis.era5Pressures(); np = length(pre); ndt = length(c1)
	ω = zeros(np,ndt)

	for idt = 1 : ndt, ip = 1 : np
		iptrop = ptrop[idt]
		isp = sp[idt]
		if pre[ip] >= ptrop[idt]
			ω[ip,idt] = c1[idt] * sin((pre[ip]-iptrop)*pi/(isp-iptrop)) + 
				        c2[idt] * sin((pre[ip]-iptrop)*2pi/(isp-iptrop))
		end
	end

	return ω,pre

end

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Reconstruction of $\omega$ from Vertical Mode Coefficients
"

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
	ωds = read_climatology(armsite,e5ds,eω)
	dt  = ωds["valid_time"][:]
	lvl = ωds["pressures"][:]
	ω   = ωds[eω.ncID][:,:]
	close(ωds)
	md"Loading Pressure Velocity data"
end

# ╔═╡ 2b99c24a-75bc-4f5b-9171-a6d297d74722
begin
	dsc = read_climatology(armsite,e5ds,eω1)
	ω1  = dsc[eω1.ncID][:]
	close(dsc)
	dsc = read_climatology(armsite,e5ds,eω2)
	ω2  = dsc[eω2.ncID][:]
	close(dsc)
	md"Loading 1st and 2nd Vertical Mode Coefficients of ω"
end

# ╔═╡ 483383ec-7356-4ad0-a5fb-abb8bf982155
begin
	dsp = read_climatology(armsite,e5ds,esp)
	sp  = dsp[esp.ncID][:] ./ 100
	close(dsp)
	dsp = read_climatology(armsite,e5ds,eptrop)
	ptrop = dsp[eptrop.ncID][:] ./ 100
	close(dsp)
	md"Loading Surface and Tropopause Pressures"
end

# ╔═╡ d3ea47b6-0570-4ad5-b3a2-c48a338c5dc6
begin
	rω,rp = ωreconstruction(ω1,ω2,ptrop,sp);
	md"Reconstructing ω from 1st and 2nd vertical modes of pressure velocity"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2001,6,3)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2001,6,8)

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
	ii = (dt .>= dtbeg) .& (dt .<= dtend)
	f1 = Figure()
	
	ax1_1 = Axis(
		f1[1,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	c1 = heatmap!(ax1_1,1:sum(ii),lvl,ω[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	lines!(ax1_1,1:sum(ii),ptrop[ii])
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax1_1,1000,10)
	
	ax1_2 = Axis(
		f1[2,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax1_2,1:sum(ii),rp,rω[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	lines!(ax1_2,1:sum(ii),ptrop[ii])
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax1_2,1000,10)
	
	ax1_3 = Axis(
		f1[3,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax1_3,1:sum(ii),rp,(rω[:,ii].-ω[:,ii])',colorrange=(-2,2),colormap=:RdBu)
	lines!(ax1_3,1:sum(ii),ptrop[ii])
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax1_3,1000,10)

	Colorbar(f1[1:3,2], c1,)
	resize_to_layout!(f1)
	f1
end

# ╔═╡ 52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
md"
### D. What Happens if the Data is Smoothed?
"

# ╔═╡ 8e2035d1-b0a1-44a0-a561-7080d5de1d09
days = 1

# ╔═╡ 55359eaf-c193-483b-8f59-6ba83a0e636a
begin
	sωds = read_climatology(armsite,e5ds,eω,days=days)
	sdt  = sωds["valid_time"][:]
	sω   = sωds[eω.ncID][:,:]; iNaN = .!isnan.(sω[1,:])
	sdt  = sdt[iNaN]
	sω   = sω[:,iNaN]
	close(sωds)
	md"Loading Smoothed Pressure Velocity data"
end

# ╔═╡ 92ecafad-296e-4441-9ef3-451093c2c560
begin
	sdsc = read_climatology(armsite,e5ds,eω1,days=days)
	sω1  = sdsc[eω1.ncID][iNaN]
	close(sdsc)
	sdsc = read_climatology(armsite,e5ds,eω2,days=days)
	sω2  = sdsc[eω2.ncID][iNaN]
	close(sdsc)
	md"Loading Smoothed 1st and 2nd Vertical Mode Coefficients of ω"
end

# ╔═╡ fa6b3cbd-f3d1-47bb-9871-158d839ae611
begin
	sdsp = read_climatology(armsite,e5ds,esp,days=days)
	ssp  = sdsp[esp.ncID][iNaN] ./ 100
	close(sdsp)
	sdsp = read_climatology(armsite,e5ds,eptrop,days=days)
	sptrop = sdsp[eptrop.ncID][iNaN] ./ 100
	close(sdsp)
	md"Loading Smoothed Troposphere Pressure Data ..."
end

# ╔═╡ 3dbd06a9-8aaa-4010-8a0c-85f169bbb8fb
begin
	srω,srp = ωreconstruction(sω1,sω2,sptrop,ssp);
	md"Reconstructing smoothed ω from 1st and 2nd vertical modes of pressure velocity"
end

# ╔═╡ 734f13a6-9c0a-4da3-8b37-4345facd7c4c
begin
	jj = (sdt .>= dtbeg) .& (sdt .<= dtend)
	f2 = Figure()
	
	ax2_1 = Axis(
		f2[1,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		yscale=log10,yticks=[10,100,1000]
	)

	c2 = heatmap!(ax2_1,1:sum(jj),lvl,sω[:,jj]',colorrange=(-0.5,0.5),colormap=:RdBu)
	lines!(ax2_1,1:sum(jj),sptrop[jj])
	ylims!(ax2_1,1000,10)
	
	ax2_2 = Axis(
		f2[2,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax2_2,1:sum(jj),rp,srω[:,jj]',colorrange=(-0.5,0.5),colormap=:RdBu)
	lines!(ax2_2,1:sum(jj),sptrop[jj])
	ylims!(ax2_2,1000,10)
	
	ax2_3 = Axis(
		f2[3,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax2_3,1:sum(jj),rp,(srω[:,jj].-sω[:,jj])',colorrange=(-0.5,0.5),colormap=:RdBu)
	lines!(ax2_3,1:sum(jj),sptrop[jj])
	ylims!(ax2_3,1000,10)

	Colorbar(f2[1:3,2], c2,)
	resize_to_layout!(f2)
	f2
end

# ╔═╡ f1208a2f-47ee-42be-9791-0389a22b9b0d
md"
### E. Decomposition of $\omega$ into Principal Components
"

# ╔═╡ 0ac5abd0-00cf-459b-abac-673717cd1b73
M = fit(PCA, ω; maxoutdim=2);

# ╔═╡ c60abbd1-f294-4500-b147-38b0ebf0e0fd
evec = eigvecs(M);

# ╔═╡ 754b0c91-4b8e-44b4-a177-23dcef174234
Y = predict(M,ω);

# ╔═╡ a7881efe-d51c-4862-b127-b47521008958
rω_pca = reconstruct(M,Y);

# ╔═╡ 15d830de-0fc7-4d09-8cdd-0081d8b17e0b
begin
	f3 = Figure()
	
	ax3_1 = Axis(
		f3[1,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	c3 = heatmap!(ax3_1,1:sum(ii),lvl,ω[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax3_1,1000,10)
	
	ax3_2 = Axis(
		f3[2,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax3_2,1:sum(ii),lvl,rω_pca[:,ii]',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax3_2,1000,10)
	
	ax3_3 = Axis(
		f3[3,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax3_3,1:sum(ii),lvl,(rω_pca[:,ii].-ω[:,ii])',colorrange=(-2,2),colormap=:RdBu)
	# xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax3_3,1000,10)

	ax3_4 = Axis(
		f3[2,2],width=50,height=200,
		# yscale=log10,yticks=[0.1]
	)

	lines!(ax3_4,evec[:,1]*0.5,lvl)
	lines!(ax3_4,evec[:,2]*0.25,lvl)
	xlims!(ax3_4,-0.2,0.2)
	ylims!(ax3_4,1000,10)

	ax3_5 = Axis(
		f3[3,2],width=200,height=200,
		# yscale=log10,yticks=[0.1]
	)

	bin1 = -15:25
	bin2 = -25:15
	h = fit(Histogram,(Y[1,:],Y[2,:]),(bin1,bin2))
	fmat = Float64.(h.weights)
	fmat[iszero.(fmat)] .= NaN
	heatmap!(ax3_5,bin1,bin2,log10.(fmat))
	# scatter!(ax3_5,Y[1,:],Y[2,:])

	# Colorbar(f3[1:3,2], c3,)
	resize_to_layout!(f3)
	f3
end

# ╔═╡ 2dfd4e50-c9d5-40a6-be3b-d9ec2bf9c4ea
md"
### F. Principal Components of Smoothed Data
"

# ╔═╡ 21527955-d742-419d-b2e3-53a7cf909158
sM = fit(PCA, sω; maxoutdim=2);

# ╔═╡ 61f401ef-182e-4298-a15d-61381f3b861c
sevec = eigvecs(sM);

# ╔═╡ 3c967cf1-a2bc-45d6-9cde-c051a6c8e747
sY = predict(sM,sω);

# ╔═╡ ccc2aa7b-7645-492a-be5e-f28f74a93021
srω_pca = reconstruct(sM,sY);

# ╔═╡ d5c3a228-b203-4a4c-a0ca-0d9af9df977f
begin
	f4 = Figure()
	
	ax4_1 = Axis(
		f4[1,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	c4 = heatmap!(ax4_1,1:sum(jj),lvl,sω[:,jj]',colorrange=(-1,1),colormap=:RdBu)
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax4_1,1000,10)
	
	ax4_2 = Axis(
		f4[2,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax4_2,1:sum(jj),lvl,srω_pca[:,jj]',colorrange=(-1,1),colormap=:RdBu)
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax4_2,1000,10)
	
	ax4_3 = Axis(
		f4[3,1],width=500,height=200,ylabel="Pressure / hPa",
		xticksvisible=false,xticklabelsvisible=false,
		# yscale=log10,yticks=[10,100,1000]
	)

	heatmap!(ax4_3,1:sum(jj),lvl,(srω[:,jj].-srω_pca[:,jj])',colorrange=(-1,1),colormap=:RdBu)
	# xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax4_3,1000,10)

	ax4_4 = Axis(
		f4[2,2],width=50,height=200,
		# yscale=log10,yticks=[0.1]
	)

	lines!(ax4_4,sevec[:,1]*0.5,lvl)
	lines!(ax4_4,sevec[:,2]*0.25,lvl)
	ylims!(ax4_4,1000,10)
	xlims!(ax4_4,-0.2,0.2)

	ax4_5 = Axis(
		f4[3,2],width=200,height=200,
		# yscale=log10,yticks=[0.1]
	)

	sbin1 = (-3:0.2:5)
	sbin2 = (-5:0.2:3)
	sh = fit(Histogram,(sY[1,:],sY[2,:]),(sbin1,sbin2))
	sfmat = Float64.(sh.weights)
	sfmat[iszero.(sfmat)] .= NaN
	heatmap!(ax4_5,sbin1,sbin2,log10.(sfmat))
	# scatter!(ax4_5,Y[1,:],Y[2,:])

	# Colorbar(f4[1:3,2], c4)
	resize_to_layout!(f4)
	f4
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
# ╟─73d960be-cbd9-4705-b68d-64bca44478ae
# ╟─aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╟─b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
# ╟─2a9efd0b-1e82-4c08-8239-451ab6ebcea0
# ╟─43502ac5-f6be-4858-8969-4e3b2e9498cd
# ╠═42aa31cf-796b-4b92-8eee-1e53c9d96765
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─2b99c24a-75bc-4f5b-9171-a6d297d74722
# ╟─483383ec-7356-4ad0-a5fb-abb8bf982155
# ╟─d3ea47b6-0570-4ad5-b3a2-c48a338c5dc6
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╠═9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
# ╠═8e2035d1-b0a1-44a0-a561-7080d5de1d09
# ╟─55359eaf-c193-483b-8f59-6ba83a0e636a
# ╟─92ecafad-296e-4441-9ef3-451093c2c560
# ╟─fa6b3cbd-f3d1-47bb-9871-158d839ae611
# ╟─3dbd06a9-8aaa-4010-8a0c-85f169bbb8fb
# ╟─734f13a6-9c0a-4da3-8b37-4345facd7c4c
# ╟─f1208a2f-47ee-42be-9791-0389a22b9b0d
# ╟─0ac5abd0-00cf-459b-abac-673717cd1b73
# ╟─c60abbd1-f294-4500-b147-38b0ebf0e0fd
# ╟─754b0c91-4b8e-44b4-a177-23dcef174234
# ╟─a7881efe-d51c-4862-b127-b47521008958
# ╟─15d830de-0fc7-4d09-8cdd-0081d8b17e0b
# ╟─2dfd4e50-c9d5-40a6-be3b-d9ec2bf9c4ea
# ╟─21527955-d742-419d-b2e3-53a7cf909158
# ╟─61f401ef-182e-4298-a15d-61381f3b861c
# ╟─3c967cf1-a2bc-45d6-9cde-c051a6c8e747
# ╟─ccc2aa7b-7645-492a-be5e-f28f74a93021
# ╟─d5c3a228-b203-4a4c-a0ca-0d9af9df977f
