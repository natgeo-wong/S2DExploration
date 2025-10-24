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
	using DSP
	using CairoMakie, LaTeXStrings
	set_theme!(theme_latexfonts())

	include(srcdir("common.jl"))
	md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 03d. ERA5 Reanalysis Power Spectra
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
### B. Defining and Loading ERA5 Datasets and Variables
"

# ╔═╡ 0394c200-c422-4aa1-9ab9-ab612e1a8019
@bind evartype Select([
	"esgl"  => "Single-Level Variable",
	"epre"  => "Pressure-Level Variable",
])

# ╔═╡ 24fba083-8237-4553-aadc-3768e6097e42
if evartype == "esgl"
	@bind ID Select([
		"blh"  => "(blh) Boundary Layer Height",
		"skt"  => "(skt) Skin Temperature",
		"t2m"  => "(t2m) 2-meter Temperature",
		"d2m"  => "(d2m) 2-meter Dewpoint",
		"tcc"  => "(tcc) Total Cloud Cover",
		"hcc"  => "(hcc) High Cloud Cover",
		"mcc"  => "(mcc) Medium Cloud Cover",
		"lcc"  => "(lcc) Low Cloud Cover",
		"sshf" => "(sshf) Surface Surface Heat Flux",
		"slhf" => "(slhf) Surface Latent Heat Flux",
		"sp"   => "(sp) Surface Pressure",
		"tciw" => "(tciw) Total Column Water Ice",
		"tclw" => "(tclw) Total Column Water Liquid",
		"tcwv" => "(tcwv) Total Column Water Vapour",
		"tcw"  => "(tcw)  Total Column Water",
		"tp"   => "(tp) Total Precipitation",
		"u10"  => "(u10) 10-meter Zonal Wind",
		"v10"  => "(v10) 10-meter Meridional Wind",
	])
else
	@bind ID Select([
		"w"  => "(w) Pressure Velocity",
		"t"  => "(t) Air Temperature",
		"q"  => "(q) Specific Humidity",
		"d"  => "(Δ) Divergence",
		"cc"  => "(cc) Cloud Cover",
	])
end

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ 953e4383-9996-4ea4-b417-44df86645d25
function load_climatology(armsite,e5ds,evar::SingleLevel;month::Int=0)

	eds = read_climatology(armsite,e5ds,evar)
	dt  = eds["valid_time"][:]
	if evartype == "esgl"
		var = eds[evar.ncID][:]
	else
		lvl = eds["pressures"][:]
		var = eds[evar.ncID][:,:]
	end
	close(eds)
	
	return var .- mean(var)

end

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
	var  = [
		load_climatology(armsite,e5ds,SingleVariable("tcc")),
		load_climatology(armsite,e5ds,SingleVariable("tcwv")),
		load_climatology(armsite,e5ds,SingleVariable("skt")),
		load_climatology(armsite,e5ds,SingleVariable("hcc")),
		load_climatology(armsite,e5ds,SingleVariable("mcc")),
		load_climatology(armsite,e5ds,SingleVariable("lcc")),
		load_climatology(armsite,e5ds,SingleVariable("ttr")),
		load_climatology(armsite,e5ds,SingleVariable("blh")),
		load_climatology(armsite,e5ds,SingleVariable("tp"))
	]
	md"Loading Single-Level Variables"
end

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Calculating the Power Spectra
"

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
	# pdg = periodogram(var,fs=24)
	# signalpower = pdg.power[2:end]
	# signalfreq  = pdg.freq[2:end]
	# ii = signalfreq .<= 1
	# signalfreq = signalfreq[ii]
	# signalpower = signalpower[ii]
	md"Calculating the Power Spectrum ..."
end

# ╔═╡ afe3e3a3-b19a-4673-a43d-9ca239a85da5
begin
	fig = Figure()

	axs = Matrix{Axis}(undef,3,3)

	for irow = 1 : 3, icol = 1 : 3
		axs[irow,icol] = Axis(
			fig[irow,icol],width=200,height=120,
			xscale=log10,yscale=log10,
			# yticklabelsvisible=isone(icol),
			xticklabelsvisible=irow==3,
		)
		xlims!(axs[irow,icol],10^-4,10)

		pdg = periodogram(var[icol+(irow-1)*3],fs=24)
		signalpower = pdg.power[2:end]
		signalfreq  = pdg.freq[2:end]
		signalpower = signalpower ./ 10. ^mean(log10.(signalpower[signalfreq.<(1/365)]))
		ymin = 1e-5
		ymax = 1e5
		ylims!(axs[irow,icol],ymin,ymax)

		lines!(axs[irow,icol],signalfreq,signalpower)
		lines!(axs[irow,icol],ones(2)/365,[ymin,ymax],color=:gray,linestyle=:dash,linewidth=0.5)
		lines!(axs[irow,icol],ones(2),[ymin,ymax],color=:gray,linestyle=:dash,linewidth=0.5)
	
		text!(10. ^-2.6,ymin*sqrt(10),text="Annual Cycle",align=(:left,:bottom),rotation=pi/2)
		text!(10. ^-0.05,ymin*sqrt(10),text="Diurnal Cycle",align=(:left,:bottom),rotation=pi/2)
	end

	Label(fig[1,1],halign=:right,valign=:top,padding=(0,5,0,5),"Total Cloud Cover")
	Label(fig[1,2],halign=:right,valign=:top,padding=(0,5,0,5),"High Cloud Cover")
	Label(fig[1,3],halign=:right,valign=:top,padding=(0,5,0,5),"OLR")
	Label(fig[2,1],halign=:right,valign=:top,padding=(0,5,0,5),"Column Water Vapor")
	Label(fig[2,2],halign=:right,valign=:top,padding=(0,5,0,5),"Medium Cloud Cover")
	Label(fig[2,3],halign=:right,valign=:top,padding=(0,5,0,5),"Boundary Layer Height")
	Label(fig[3,1],halign=:right,valign=:top,padding=(0,5,0,5),"Skin Temperature")
	Label(fig[3,2],halign=:right,valign=:top,padding=(0,5,0,5),"Low Cloud Cover")
	Label(fig[3,3],halign=:right,valign=:top,padding=(0,5,0,5),"Rainfall Rate")

	Label(fig[0,:],"Power Spectra for ERA5 Variables at $armsite ARM Site (1980-2024)",font=:bold)
	Label(fig[4,:],L"Frequency $f$ / day$^{-1}$")
	Label(fig[:,0],L"A$^2f$",rotation=pi/2)
	
	resize_to_layout!(fig)
	fig
end

# ╔═╡ b80a1b46-ce7b-4d75-bcd3-4a1ebca4e776
CairoMakie.save(plotsdir("03-powerspectra-$armsite.png"),fig)

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─0394c200-c422-4aa1-9ab9-ab612e1a8019
# ╟─24fba083-8237-4553-aadc-3768e6097e42
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─953e4383-9996-4ea4-b417-44df86645d25
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─afe3e3a3-b19a-4673-a43d-9ca239a85da5
# ╟─b80a1b46-ce7b-4d75-bcd3-4a1ebca4e776
