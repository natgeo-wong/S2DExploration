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
# 03a. Exploring the Diurnal Cycle
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

# ╔═╡ 24fba083-8237-4553-aadc-3768e6097e42
@bind ID_sgl Select([
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

# ╔═╡ 5edcc0ff-6ddc-470a-88c9-b980545152d6
@bind ID_pre Select([
    "w"  => "(w) Pressure Velocity",
    "t"  => "(t) Air Temperature",
    "q"  => "(q) Specific Humidity",
])

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
esgl = SingleVariable(ID_sgl)

# ╔═╡ b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
epre = PressureVariable(ID_pre)

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
function load_climatology(armsite,e5ds,evar::SingleLevel;mo::Int)

    sds = read_climatology(armsite,e5ds,esgl)
    dt  = sds["valid_time"][:]
    ii = month.(dt) .== mo
    sgl = sds[esgl.ncID][ii]
    close(sds)
    sgl = dropdims(mean(reshape(sgl,24,:),dims=2),dims=2)

    slon,_ = getarmcoordinates()
    t_s,sgl = utc2local(sgl,slon)
    
    return t_s,sgl

end

# ╔═╡ 9ed3ced3-a38c-4edf-9cfd-510b564661de
function load_climatology(armsite,e5ds,evar::PressureLevel;mo::Int)

    pds = read_climatology(armsite,e5ds,epre)
    dt  = pds["valid_time"][:]
    lvl = pds["pressures"][:]; nlvl = length(lvl)
    ii = month.(dt) .== mo
    pre = pds[evar.ncID][:,ii]
    close(pds)
    pre = dropdims(mean(reshape(pre,nlvl,24,:),dims=3),dims=3)
    
    slon,_ = getarmcoordinates()
    t_p,pre = utc2local(pre,slon)
    
    return t_p,lvl,pre

end

# ╔═╡ ca6c960a-a8c6-43f0-86d8-841e0f520b33
md"Month: $(@bind mo PlutoUI.Slider(1:12,default=1, show_value=true))"

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
    t_s,sgl = load_climatology(armsite,e5ds,esgl,mo=mo)
    md"Loading Single-Level Variable: $(esgl.ID)"
end

# ╔═╡ 2b99c24a-75bc-4f5b-9171-a6d297d74722
begin
    t_p,lvl,pre = load_climatology(armsite,e5ds,epre,mo=mo)
    md"Loading Pressure-Level Variable: $(epre.ID)"
end

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
    f1 = Figure()
    
    ax1_1 = Axis(
        f1[1,1],width=500,height=300,ylabel="Pressure / hPa",
        xticklabelsvisible=false,
        # yscale=log10,yticks=[10,100,1000]
    )

    c1 = heatmap!(ax1_1,t_p,lvl,pre',colorrange=(-1,1).*1e-1,colormap=:RdBu)
    xlims!(ax1_1,0,24)
    ylims!(ax1_1,1000,10)

    # axislegend()
    
    ax1_2 = Axis(f1[2,1],width=500,height=50,xlabel="Local Time")

    lines!(ax1_2,t_s,sgl)
    xlims!(ax1_2,0,24)

    Colorbar(f1[1,2], c1,)
    resize_to_layout!(f1)
    f1
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─2503e5bc-830c-4a40-a0f8-65ec8e669bb5
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─24fba083-8237-4553-aadc-3768e6097e42
# ╟─5edcc0ff-6ddc-470a-88c9-b980545152d6
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╟─b5741a76-cbdd-4f0c-bd22-e9de5cd9ceea
# ╟─31b67b50-420e-4c8e-8cf9-1d4b0e69137e
# ╟─ddd75fee-7fbc-4756-bb98-8b1bb3a888f6
# ╟─fd353a8b-72a3-4c58-97f1-bb8e94bbbcf4
# ╟─a3a840c6-062c-48a7-a66d-1ec5e375e671
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─10662679-91f7-43e6-83fd-b9eeb02794ca
# ╟─9ed3ced3-a38c-4edf-9cfd-510b564661de
# ╟─ca6c960a-a8c6-43f0-86d8-841e0f520b33
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─2b99c24a-75bc-4f5b-9171-a6d297d74722
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
