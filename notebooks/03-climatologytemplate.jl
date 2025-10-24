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
	using DataInterpolations
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())

	include(srcdir("smoothing.jl"))
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 03. Climatology Figure Template
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

# ╔═╡ 2503e5bc-830c-4a40-a0f8-65ec8e669bb5
getarmcoordinates() = if armsite == "SGP"
	return -97.487643, 36.607322;
else
	return -87.338177, 34.342481;
end

# ╔═╡ e282e373-f8a9-4a11-b1a7-e1d270b14090
md"
### B. Defining ERA5 Datasets and Variables
"

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ 31b67b50-420e-4c8e-8cf9-1d4b0e69137e
md"
### C. Defining Interpolation Functions to Local Time
"

# ╔═╡ ddd75fee-7fbc-4756-bb98-8b1bb3a888f6
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ fd353a8b-72a3-4c58-97f1-bb8e94bbbcf4
function utc2local(data :: Vector,longitude)

	nt = length(data)
	t = 0:23
	ts = longitude2timeshift.(longitude)
	it = (0:0.5:24); it = (it[2:end].+it[1:(end-1)])/2; nit = length(it)

	tl = t .+ ts

	itp = AkimaInterpolation(vcat(data,data,data), vcat(tl.-24,tl,tl.+24))
	var = itp.(it)

	return it,var

end

# ╔═╡ a3a840c6-062c-48a7-a66d-1ec5e375e671
function utc2local(data :: Matrix,longitude)

	np,nt = size(data)
	t = 0:23
	ts = longitude2timeshift.(longitude)
	it = (0:0.5:24); it = (it[2:end].+it[1:(end-1)])/2; nit = length(it)

	var = zeros(np,nit)

	tl = t .+ ts

	for ip = 1 : np
		idata = @views data[ip,:]
		itp = AkimaInterpolation(vcat(idata,idata,idata), vcat(tl.-24,tl,tl.+24))
		var[ip,:] = itp.(it)
	end

	return it,var

end

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### D. Loading Timeseries Data for Stations
"

# ╔═╡ 10662679-91f7-43e6-83fd-b9eeb02794ca
function load_climatology(armsite,e5ds,evar::SingleLevel;month::Int=0)

	sds = read_climatology(armsite,e5ds,evar)
	dt  = sds["valid_time"][:]
	jj = dt .>= Date(1990)
	if iszero(month)
		sgl = sds[evar.ncID][jj]
	else
		ii = Month(dt) .== month
		sgl = sds[evar.ncID][jj][ii]
	end
	close(sds)
	sgl = dropdims(mean(reshape(sgl,24,:),dims=2),dims=2)

	slon,_ = getarmcoordinates()
	t_s,sgl = utc2local(sgl,slon)
	
	return t_s,sgl

end

# ╔═╡ 9ed3ced3-a38c-4edf-9cfd-510b564661de
function load_climatology(armsite,e5ds,evar::PressureLevel;month::Int=0)

	pds = read_climatology(armsite,e5ds,evar)
	dt  = pds["valid_time"][:]
	jj = dt .>= Date(1990)
	lvl = pds["pressures"][:]; nlvl = length(lvl)
	if iszero(month)
		pre = pds[evar.ncID][:,jj]
	else
		ii = Month(dt) .== month
		pre = pds[evar.ncID][:,jj][:,ii]
	end
	close(pds)
	pre = dropdims(mean(reshape(pre,nlvl,24,:),dims=3),dims=3)
	
	slon,_ = getarmcoordinates()
	t_p,pre = utc2local(pre,slon)
	
	return t_p,lvl,pre

end

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
	t_s,tcc = load_climatology(armsite,e5ds,SingleVariable("tcc"))
	_,tcwv  = load_climatology(armsite,e5ds,SingleVariable("tcwv"))
	_,skt   = load_climatology(armsite,e5ds,SingleVariable("skt"))
	_,hcc   = load_climatology(armsite,e5ds,SingleVariable("hcc"))
	_,mcc   = load_climatology(armsite,e5ds,SingleVariable("mcc"))
	_,lcc   = load_climatology(armsite,e5ds,SingleVariable("lcc"))
	_,ttr   = load_climatology(armsite,e5ds,SingleVariable("ttr"))
	_,blh   = load_climatology(armsite,e5ds,SingleVariable("blh"))
	_,tp    = load_climatology(armsite,e5ds,SingleVariable("tp"))
	md"Loading single level variable climatologies ..."
end

# ╔═╡ 2b99c24a-75bc-4f5b-9171-a6d297d74722
begin
	t_p,lvl,w = load_climatology(armsite,e5ds,PressureVariable("w"));
	_,_,t = load_climatology(armsite,e5ds,PressureVariable("t"));
	_,_,q = load_climatology(armsite,e5ds,PressureVariable("q"));
	md"Loading pressure level variable climatologies ..."
end

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
	fig = Figure()

	axs = Matrix{Axis}(undef,4,3)

	for irow = 1 : 4, icol = 1 : 3
		if isone(irow)
			axs[irow,icol] = Axis(
				fig[irow,icol],width=200,height=120,
				xticklabelsvisible=false,
				# yticklabelsvisible=isone(icol),
				# yscale=log10,yticks=[10,100,1000]
			)
			ylims!(axs[irow,icol],1000,10)
		else
			axs[irow,icol] = Axis(
				fig[irow+1,icol],width=200,height=50,
				xticklabelsvisible=irow==4,
				# yticklabelsvisible=isone(icol),
			)
		end
		xlims!(axs[irow,icol],0,24)
	end
	

	cbr_1 = heatmap!(axs[1,1],t_p,lvl,w',colorrange=(-1,1).*1e-1,colormap=:RdBu)
	cbr_2 = heatmap!(axs[1,2],t_p,lvl,(t.-mean(t,dims=2))',colorrange=(-2,2),colormap=:balance)
	cbr_3 = heatmap!(axs[1,3],t_p,lvl,(q.-mean(q,dims=2))'*1000,colorrange=(-5,5).*1e-1,colormap=:BrBg)

	lines!(axs[2,1],t_s,tcc*100)
	lines!(axs[3,1],t_s,tcwv)
	lines!(axs[4,1],t_s,skt)
	lines!(axs[2,2],t_s,hcc*100)
	lines!(axs[3,2],t_s,mcc*100)
	lines!(axs[4,2],t_s,lcc*100)
	lines!(axs[2,3],t_s,-ttr/3600)
	lines!(axs[3,3],t_s,blh/1000)
	lines!(axs[4,3],t_s,tp*1000)

	axs[1,1].ylabel = "Pressure / hPa"
	axs[2,1].ylabel = "%"
	axs[3,1].ylabel = "mm"
	axs[4,1].ylabel = "K"
	axs[2,2].ylabel = "%"
	axs[3,2].ylabel = "%"
	axs[4,2].ylabel = "%"
	axs[2,3].ylabel = L"W m$^{-2}$"
	axs[3,3].ylabel = "km"
	axs[4,3].ylabel = L"mm hr$^{-1}$"

	Label(fig[1,1],halign=:right,valign=:top,padding=(0,5,0,5),L"\omega")
	Label(fig[1,2],halign=:right,valign=:top,padding=(0,7,0,5),L"T'")
	Label(fig[1,3],halign=:right,valign=:top,padding=(0,7,0,5),L"q'")
	Label(fig[3,1],halign=:right,valign=:top,padding=(0,5,0,5),"Total Cloud Cover")
	Label(fig[3,2],halign=:right,valign=:top,padding=(0,5,0,5),"High Cloud Cover")
	Label(fig[3,3],halign=:right,valign=:top,padding=(0,5,0,5),"OLR")
	Label(fig[4,1],halign=:right,valign=:top,padding=(0,5,0,5),"Column Water Vapor")
	Label(fig[4,2],halign=:right,valign=:top,padding=(0,5,0,5),"Medium Cloud Cover")
	Label(fig[4,3],halign=:right,valign=:top,padding=(0,5,0,5),"Boundary Layer Height")
	Label(fig[5,1],halign=:right,valign=:top,padding=(0,5,0,5),"Skin Temperature")
	Label(fig[5,2],halign=:right,valign=:top,padding=(0,5,0,5),"Low Cloud Cover")
	Label(fig[5,3],halign=:right,valign=:top,padding=(0,5,0,5),"Rainfall Rate")

	Label(fig[0,:],"Climatology at the $armsite ARM Site (1980-2024)",font=:bold)
	Label(fig[6,:],"Hour of Day")

	Colorbar(fig[2,1],cbr_1,vertical=false,flipaxis=false,label=L"Pa s$^{-1}$")
	Colorbar(fig[2,2],cbr_2,vertical=false,flipaxis=false,label="K")
	Colorbar(fig[2,3],cbr_3,vertical=false,flipaxis=false,label=L"g kg$^{-1}$")
	resize_to_layout!(fig)
	fig
end

# ╔═╡ a6d6a109-fe25-417d-8334-2df5cacfe318
CairoMakie.save(plotsdir("03-climatology-$armsite.png"),fig)

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─2503e5bc-830c-4a40-a0f8-65ec8e669bb5
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─31b67b50-420e-4c8e-8cf9-1d4b0e69137e
# ╟─ddd75fee-7fbc-4756-bb98-8b1bb3a888f6
# ╟─fd353a8b-72a3-4c58-97f1-bb8e94bbbcf4
# ╟─a3a840c6-062c-48a7-a66d-1ec5e375e671
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─10662679-91f7-43e6-83fd-b9eeb02794ca
# ╟─9ed3ced3-a38c-4edf-9cfd-510b564661de
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─2b99c24a-75bc-4f5b-9171-a6d297d74722
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─a6d6a109-fe25-417d-8334-2df5cacfe318
