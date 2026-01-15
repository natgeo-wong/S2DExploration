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
    using NASAPrecipitation
    using CairoMakie, LaTeXStrings
    set_theme!(theme_latexfonts())

    include(srcdir("smoothing.jl"))
    md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02b. Exploration of ERA5 Reanalysis Timeseries Data
"

# ╔═╡ 370903b1-1c85-4337-89c0-fd7cc358359d
TableOfContents()

# ╔═╡ 2205ab68-6602-4222-b3d5-1ebf2a7b8082
md"
### A. Loading Station Information
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

# ╔═╡ e282e373-f8a9-4a11-b1a7-e1d270b14090
md"
### B. Defining IMERG and ERA5 Datasets and Variables
"

# ╔═╡ 326187cb-6cf9-49bb-99f7-a27d1c825eae
npd = IMERGFinalHH(start=Date(2001),stop=Date(2024,12),path=datadir())

# ╔═╡ c3329018-a49d-4bc3-b6d3-5bc4ea401157
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
esgl = SingleVariable("tp")

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Timeseries Data for Stations
"

# ╔═╡ 6cbdd12a-1cd4-44b4-8793-1e95ed5c984e
begin
    ids = read_climatology(armsite,npd)
    npd_dt = ids["time"][:]
    npd_tp = ids["precipitation"][:]
    close(ids)
    md"Loading NASA Precipitation Dataset"
end

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
    sds = read_climatology(armsite,e5ds,esgl)
    era_dt = sds["valid_time"][:]
    era_tp = sds[esgl.ncID][:]
    close(sds)
    md"Loading Single-Level Variable: $(esgl.ID)"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2001,6)

# ╔═╡ c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
dtend = Date(2001,6,10)

# ╔═╡ 52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
md"
### D. What Happens if the Data is Smoothed?
"

# ╔═╡ 8e2035d1-b0a1-44a0-a561-7080d5de1d09
days = 30

# ╔═╡ 92ecafad-296e-4441-9ef3-451093c2c560
begin
    sera_tp = deepcopy(era_tp); smoothing!(sera_tp,days=days)
    md"Smoothing Single-Level Data ..."
end

# ╔═╡ fa6b3cbd-f3d1-47bb-9871-158d839ae611
begin
    snpd_tp = deepcopy(npd_tp); smoothing!(snpd_tp,days=days*2)
    md"Smoothing IMERG Precipitation data ..."
end

# ╔═╡ 6c3b781d-e602-437e-b312-8f9d48bd1cc5
md"
### E. Raw minus Smoothed
"

# ╔═╡ 3159caf0-0443-4a64-8001-82799e3c4865
# ╠═╡ disabled = true
#=╠═╡
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

    # axislegend()
    
    ax3_2 = Axis(f3[2,1],width=500,height=50,xlabel="Date")

    lines!(ax3_2,dt[ii],sgl[ii])
    lines!(ax3_2,dt[ii],ssgl[ii])
    xlims!(ax3_2,dtbeg,dtend)

    Colorbar(f3[1,2], c3,)
    resize_to_layout!(f3)
    f3
end
  ╠═╡ =#

# ╔═╡ 734f13a6-9c0a-4da3-8b37-4345facd7c4c
begin
    f2 = Figure()
    
    ax2_1 = Axis(
        f2[1,1],width=500,height=50,ylabel=xlabel="Date"
    )

    lines!(ax2_1,npd_dt[ii],npd_tp[ii])
    lines!(ax2_1,era_dt[ii],era_tp[ii])
    xlims!(ax2_1,dtbeg,dtend)

    # axislegend()

    resize_to_layout!(f2)
    f2
end

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
begin
    ii = (dt .>= dtbeg) .& (dt .<= dtend)
    f1 = Figure()
    
    ax1_1 = Axis(
        f1[1,1],width=500,height=50,ylabel=xlabel="Date"
    )

    lines!(ax1_1,npd_dt[ii],npd_tp[ii])
    lines!(ax1_1,era_dt[ii],era_tp[ii])
    xlims!(ax1_1,dtbeg,dtend)

    # axislegend()

    resize_to_layout!(f1)
    f1
end

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─370903b1-1c85-4337-89c0-fd7cc358359d
# ╟─2205ab68-6602-4222-b3d5-1ebf2a7b8082
# ╟─baaf05e1-f7a8-49f1-850d-d42cb2ed0f2e
# ╟─8a43991e-cf8e-4304-9d0e-f235a0da1f36
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╠═326187cb-6cf9-49bb-99f7-a27d1c825eae
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╟─aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─6cbdd12a-1cd4-44b4-8793-1e95ed5c984e
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═c2176f38-e2e2-4bfb-85e6-5e2f77fd1cf4
# ╟─9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╟─52d92463-ad3c-42e1-a5bf-4d3ab32a2f82
# ╠═8e2035d1-b0a1-44a0-a561-7080d5de1d09
# ╟─92ecafad-296e-4441-9ef3-451093c2c560
# ╟─fa6b3cbd-f3d1-47bb-9871-158d839ae611
# ╠═734f13a6-9c0a-4da3-8b37-4345facd7c4c
# ╟─6c3b781d-e602-437e-b312-8f9d48bd1cc5
# ╟─3159caf0-0443-4a64-8001-82799e3c4865
