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
	using RegionGrids, ERA5Reanalysis, NASAMergedTb
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())

	include(srcdir("common.jl"))
	include(srcdir("calculate.jl"))
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02g. NASA vs ERA5 Comparison: Outgoing Longwave Radiation
"

# ╔═╡ 370903b1-1c85-4337-89c0-fd7cc358359d
TableOfContents()

# ╔═╡ 2205ab68-6602-4222-b3d5-1ebf2a7b8082
md"
### A. Loading Station and ERA5 Dataset Information
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

# ╔═╡ 27c17138-69a7-4fdd-bd66-8904a007d509
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ 28d526bc-bf12-48a0-82a8-bcefc38c35ae
evar = SingleVariable("ttr")

# ╔═╡ 1f7f1d77-85cd-4bac-921c-6e43d7f00991
elon,elat = ERA5Reanalysis.nativelonlat()

# ╔═╡ 149ad521-c0dc-4c09-b1d8-127c2500c23d
begin
	blon,blat = NASAMergedTb.btdlonlat()
	ggrd1 = RegionGrid(GeoRegion("SGP",path=srcdir()),blon,blat)
	nlon_SGP = length(ggrd1.lon); nlat_SGP = length(ggrd1.lat)
	blongrid1 = zeros(nlon_SGP,nlat_SGP); blongrid1 .= ggrd1.lon
	blatgrid1 = zeros(nlon_SGP,nlat_SGP); blatgrid1 .= reshape(ggrd1.lat,1,:)
	ggrd2 = RegionGrid(GeoRegion("BNF",path=srcdir()),blon,blat)
	nlon_BNF = length(ggrd2.lon); nlat_BNF = length(ggrd2.lat)
	blongrid2 = zeros(nlon_BNF,nlat_BNF); blongrid2 .= ggrd2.lon
	blatgrid2 = zeros(nlon_BNF,nlat_BNF); blatgrid2 .= reshape(ggrd2.lat,1,:);
	md"Loading Brightness Temperature Grid Information ..."
end

# ╔═╡ bba26d40-6808-485e-9834-82ea4dc6c176
md"
### B. Loading OLR Timeseries Data for merged-IR and ERA5
"

# ╔═╡ 1a48c947-198a-4e85-af2d-87dc2eefe94a
if armsite == "SGP"
	@bind iistn Select([ii => "($(sgp_info[ii,1])) $(sgp_info[ii,2])" for ii in 1 : nsgp])
else
	@bind iistn Select([ii => "($(bnf_info[ii,1])) $(bnf_info[ii,2])" for ii in 1 : nbnf])
end

# ╔═╡ 17897c14-ebce-4947-8c85-668d4818423f
if armsite == "SGP"
	stnID = sgp_info[iistn,1];
else
	stnID = bnf_info[iistn,1];
end;

# ╔═╡ d9818f36-6dd5-45d1-b56b-42392ed79049
begin
	tds = NCDataset(datadir(
		"mergedIR","timeseries-OLR",
		"OLR-$(armsite)_$(stnID)-20010101-20241231.nc"
	))
	ndt = tds["valid_time"][:]
	olr = tds["OLR"][:]
	close(tds)
end

# ╔═╡ 215d963f-3019-46ad-b139-569df1cad1e9
begin
	eds = read_climatology(armsite,e5ds,evar)
	edt = eds["valid_time"][:]
	var = eds[evar.ncID][:] / -3600
	close(eds)
	md"Loading Single-Level Variable: $(evar.ID)"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2001,7)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2001,7,10)

# ╔═╡ 634c8c07-1aac-4763-86ae-b2b13519164a
begin
	ii = (ndt .>= dtbeg) .& (ndt .<= dtend)
	jj = (edt .>= dtbeg) .& (edt .<= dtend)
	t1 = (ndt .>= Date(2001)) .& (ndt .<= Date(2024,12,31))
	t2 = (edt .>= Date(2001)) .& (edt .<= Date(2024,12,31))
	f1 = Figure()

	if armsite == "SGP"
		stnstr = sgp_info[iistn,1]
	else
		stnstr = bnf_info[iistn,1]
	end
	
	ax1_1 = Axis(
		f1[1,1],width=400,height=200,yticks=100:50:350,xlabel="Date",
		ylabel=L"OLR / W m$^{-2}$",yminorticks=50:10:350,yminorticksvisible=true,
	)

	lines!(ax1_1,ndt[ii],olr[ii],label=stnstr,color=:black,linewidth=3)
	lines!(ax1_1,edt[jj],var[jj],label=stnstr,color=:blue,linewidth=2)
	xlims!(ax1_1,dtbeg,dtend)
	ylims!(ax1_1,75,350)

	ax1 = Axis(
		f1[1,2],width=200,height=200,limits=(75,350,75,350),
		xticks=100:50:350,xlabel=L"ERA5 OLR / W m$^{-2}$",
		xminorticks=50:10:350,xminorticksvisible=true,
		yticks=100:50:350,ylabel=L"NASA $T_b$-Derived OLR / W m$^{-2}$",
		yminorticks=50:10:350,yminorticksvisible=true,yaxisposition=:right
	)

	varii = (var[t2][1:(end-1)] .+ var[t2][2:end])/2
	olrii = dropdims(mean(reshape(olr[t1],2,:),dims=1),dims=1)
	bins = 75:5:350
	h = fit(Histogram,(varii,olrii),(bins,bins))
	fmat = Float64.(h.weights)
	fmat[iszero.(fmat)] .= NaN

	c1 = heatmap!(ax1,bins,bins,fmat,)
	lines!(ax1,[50,350],[50,350],color=:grey,linewidth=0.5)

	resize_to_layout!(f1)
	f1
end

# ╔═╡ 77709edd-a599-43ab-8e38-f1d4c4b496bc
md"
### C. Intragrid Variability of OLR
"

# ╔═╡ c93925c4-8d9c-4313-bfde-3b45bed721ae
begin
	iSGP = closestnativelonlat(Point2(sgp_info[1,5],sgp_info[1,6]))
	iBNF = closestnativelonlat(Point2(bnf_info[1,5],bnf_info[1,6]))

	iarg_SGP = closestnativelonlat(Point2.(blongrid1[:],blatgrid1[:])) .== iSGP
	iarg_BNF = closestnativelonlat(Point2.(blongrid2[:],blatgrid2[:])) .== iBNF

md"
1. Find ielon/ielat closest to stations
2. Find blon/blats that are closest to ielon/ielat
"
end

# ╔═╡ 8bde2ae9-60ae-4706-b585-26a4a5c3c835
begin
	fig2 = Figure()

	ax2_1 = Axis(
			fig2[1,1],width=200,height=200,title="Southern Great Plains (New)",
			limits=(-97.7,-96.9,36.3,37.1),xticks=-100:0.5:95,yticks=34:0.5:39
	)
	
	scatter!(ax2_1,elon.-360,elat,markersize=10,color=:black)
	scatter!(ax2_1,blongrid1[:],blatgrid1[:],markersize=5,color=:grey)
	scatter!(ax2_1,blongrid1[iarg_SGP],blatgrid1[iarg_SGP],markersize=5,color=:lightblue)
	scatter!(ax2_1,elon[iSGP].-360,elat[iSGP],markersize=10,color=:blue)
	scatter!(ax2_1,sgp_info[1,5],sgp_info[1,6],markersize=15,color=:red)

	ax2_2 = Axis(
			fig2[1,2],width=200,height=200,title="Bankhead National Forest",
			limits=(-87.7,-86.9,34.0,34.8),xticks=-90:0.5:-85,yticks=32:0.5:37,
	)
	
	scatter!(ax2_2,elon.-360,elat,markersize=10,color=:black)
	scatter!(ax2_2,blongrid2[:],blatgrid2[:],markersize=5,color=:grey)
	scatter!(ax2_2,blongrid2[iarg_BNF],blatgrid2[iarg_BNF],markersize=5,color=:lightblue)
	scatter!(ax2_2,elon[iBNF].-360,elat[iBNF],markersize=10,color=:blue)
	scatter!(ax2_2,bnf_info[1,5],bnf_info[1,6],markersize=15,color=:red)

	Label(fig2[2,:],L"Longitude / $\degree$")
	Label(fig2[1,0],L"Latitude / $\degree$",rotation=pi/2)

	# Colorbar(fig2[1,4], c2, label="Topographic Height / km")
	resize_to_layout!(fig2)
	fig2

end

# ╔═╡ 5bb4c5da-1fe4-4ce4-8736-451117dc13f6
begin
	nds = NCDataset(datadir(
		"mergedIR","timeseries-Tb_compare",
		"Tb_compare-$(armsite)_$(stnID)-20010101-20241231.nc"
	))
	atb = nds["Tb"][:,:]; npnt = size(atb,1)
	μTb = dropdims(mean(atb,dims=1),dims=1)
	aolr = Tb2OLR.(atb)
	μolr = Tb2OLR.(μTb)
	close(nds)
end

# ╔═╡ 0b11ca14-4b4d-464d-8c77-c9823e9b376c
begin
	f3 = Figure()
	
	ax3_1 = Axis(
		f3[1,1],width=400,height=200,yticks=100:50:350,xlabel="Date",
		ylabel=L"OLR / W m$^{-2}$",yminorticks=50:10:350,yminorticksvisible=true,
	)

	for ipnt = 1 : npnt
		lines!(ax3_1,ndt[ii],aolr[ipnt,ii],color=:grey,linewidth=1)
	end
	lines!(ax3_1,ndt[ii],μolr[ii],label="NASA",color=:black,linewidth=3)
	lines!(ax3_1,edt[jj],var[jj],label="ERA5",color=:blue,linewidth=3)
	xlims!(ax3_1,dtbeg,dtend)
	ylims!(ax3_1,75,350)

	axislegend(position=:lb,framevisible=:false)

	ax3_2 = Axis(
		f3[1,2],width=200,height=200,limits=(75,350,75,350),
		xticks=100:50:350,xlabel=L"ERA5 OLR / W m$^{-2}$",
		xminorticks=50:10:350,xminorticksvisible=true,
		yticks=100:50:350,ylabel=L"NASA $T_b$-Derived OLR / W m$^{-2}$",
		yminorticks=50:10:350,yminorticksvisible=true,yaxisposition=:right
	)

	μolrii = dropdims(mean(reshape(μolr[t1],2,:),dims=1),dims=1)
	μh = fit(Histogram,(varii,μolrii),(bins,bins))
	μfmat = Float64.(μh.weights)
	μfmat[iszero.(μfmat)] .= NaN

	c3 = heatmap!(ax3_2,bins,bins,μfmat,)
	lines!(ax3_2,[50,350],[50,350],color=:grey,linewidth=0.5)

	resize_to_layout!(f3)
	f3
end

# ╔═╡ 5421f7e1-83a4-4e1d-b237-3f969f60e1a5
CairoMakie.save(plotsdir("02g-OLRcomparison.png"),f3)

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─baaf05e1-f7a8-49f1-850d-d42cb2ed0f2e
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─27c17138-69a7-4fdd-bd66-8904a007d509
# ╟─28d526bc-bf12-48a0-82a8-bcefc38c35ae
# ╟─1f7f1d77-85cd-4bac-921c-6e43d7f00991
# ╟─149ad521-c0dc-4c09-b1d8-127c2500c23d
# ╟─bba26d40-6808-485e-9834-82ea4dc6c176
# ╟─1a48c947-198a-4e85-af2d-87dc2eefe94a
# ╟─17897c14-ebce-4947-8c85-668d4818423f
# ╟─d9818f36-6dd5-45d1-b56b-42392ed79049
# ╟─215d963f-3019-46ad-b139-569df1cad1e9
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╟─634c8c07-1aac-4763-86ae-b2b13519164a
# ╟─77709edd-a599-43ab-8e38-f1d4c4b496bc
# ╟─c93925c4-8d9c-4313-bfde-3b45bed721ae
# ╟─8bde2ae9-60ae-4706-b585-26a4a5c3c835
# ╟─5bb4c5da-1fe4-4ce4-8736-451117dc13f6
# ╟─0b11ca14-4b4d-464d-8c77-c9823e9b376c
# ╟─5421f7e1-83a4-4e1d-b237-3f969f60e1a5
