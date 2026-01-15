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
    using ERA5Reanalysis
    using CairoMakie, LaTeXStrings
    set_theme!(theme_latexfonts())

    include(srcdir("smoothing.jl"))
    md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02d. Calculating the Tropopause Pressure Height
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

# ╔═╡ 5edcc0ff-6ddc-470a-88c9-b980545152d6
@bind ID_pre Select([
    "w"  => "(w) Pressure Velocity",
    "t"  => "(t) Air Temperature",
    "q"  => "(q) Specific Humidity",
    "dtdp"  => "(dtdp) Temperature Gradient over Pressure",
    "dtdp_hres"  => "(dtdp) High-Resolution Temperature Gradient over Pressure",
])

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
epre = PressureVariable(ID_pre,path=srcdir())

# ╔═╡ 323a3977-466e-4609-8a3c-42c0bbab230f
etrp = SingleVariable("ptrop",path=srcdir())

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Timeseries Data for Stations
"

# ╔═╡ 2b99c24a-75bc-4f5b-9171-a6d297d74722
begin
    pds = read_climatology(armsite,e5ds,epre)
    dt  = pds["valid_time"][:]
    lvl = pds["pressures"][:]
    pre = pds[epre.ncID][:,:] * 1000
    close(pds)
    md"Loading Pressure-Level Variable: $(epre.ID)"
end

# ╔═╡ 1dd94d0c-5800-47dd-a571-bcc2652187c1
begin
    tds = read_climatology(armsite,e5ds,etrp)
    ptp = tds[etrp.ncID][:] / 100
    close(tds)
    md"Loading Single-Level Variable: $(etrp.ID)"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2011)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2021)

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
    ii = (dt .>= dtbeg) .& (dt .<= dtend)
    f1 = Figure()
    
    ax1_1 = Axis(
        f1[1,1],width=500,height=300,ylabel="Pressure / hPa",
        xticksvisible=false,xticklabelsvisible=false,
        # yscale=log10,yticks=[10,100,1000]
    )

    c1 = heatmap!(ax1_1,1:sum(ii),lvl,pre[:,ii]',colorrange=(-2,2),colormap=:RdBu)
    lines!(ax1_1,1:sum(ii),ptp[ii])
    # xlims!(ax1_1,dtbeg,dtend)
    ylims!(ax1_1,1000,10)

    Colorbar(f1[1,2], c1,)
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
    spds = read_climatology(armsite,e5ds,epre,days=days)
    spre = spds[epre.ncID][:,:] * 1000
    close(spds)
    md"Loading Smoothed Troposphere Pressure Data ..."
end

# ╔═╡ fa6b3cbd-f3d1-47bb-9871-158d839ae611
begin
    stds = read_climatology(armsite,e5ds,etrp,days=days)
    ptps = stds[etrp.ncID][:] / 100
    close(stds)
    md"Loading Smoothed Troposphere Pressure Data ..."
end

# ╔═╡ 734f13a6-9c0a-4da3-8b37-4345facd7c4c
begin
    f2 = Figure()
    
    ax2_1 = Axis(
        f2[1,1],width=500,height=300,ylabel="Pressure / hPa",
        xticksvisible=false,xticklabelsvisible=false,
        # yscale=log10,yticks=[10,100,1000]
    )

    c2 = heatmap!(ax2_1,1:sum(ii),lvl,spre[:,ii]',colorrange=(-2,2),colormap=:RdBu)
    lines!(ax2_1,1:sum(ii),ptps[ii])
    # xlims!(ax1_1,dtbeg,dtend)
    ylims!(ax2_1,1000,10)

    Colorbar(f2[1,2], c2,)
    resize_to_layout!(f2)
    f2
end

# ╔═╡ 6c3b781d-e602-437e-b312-8f9d48bd1cc5
md"
### E. Raw minus Smoothed
"

# ╔═╡ 3159caf0-0443-4a64-8001-82799e3c4865
begin
    f3 = Figure()
    
    ax3_1 = Axis(
        f3[1,1],width=500,height=300,ylabel="Pressure / hPa",
        xticksvisible=false,xticklabelsvisible=false,
        # yscale=log10,yticks=[10,100,1000]
    )

    c3 = heatmap!(ax3_1,1:sum(ii),lvl,(pre[:,ii].-spre[:,ii])',colorrange=(-2,2),colormap=:RdBu)
    # xlims!(ax1_1,dtbeg,dtend)
    ylims!(ax3_1,1000,10)

    Colorbar(f3[1,2], c3,)
    resize_to_layout!(f3)
    f3
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─5edcc0ff-6ddc-470a-88c9-b980545152d6
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
# ╟─323a3977-466e-4609-8a3c-42c0bbab230f
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─2b99c24a-75bc-4f5b-9171-a6d297d74722
# ╟─1dd94d0c-5800-47dd-a571-bcc2652187c1
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╠═9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
# ╠═8e2035d1-b0a1-44a0-a561-7080d5de1d09
# ╟─55359eaf-c193-483b-8f59-6ba83a0e636a
# ╟─fa6b3cbd-f3d1-47bb-9871-158d839ae611
# ╟─734f13a6-9c0a-4da3-8b37-4345facd7c4c
# ╟─6c3b781d-e602-437e-b312-8f9d48bd1cc5
# ╟─3159caf0-0443-4a64-8001-82799e3c4865
