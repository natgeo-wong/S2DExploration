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
    using CairoMakie, LaTeXStrings
    set_theme!(theme_latexfonts())

    include(srcdir("smoothing.jl"))
    md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 04b. Top of the Boundary Layer
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

# ╔═╡ bd0593a8-c898-4707-90cc-6ddb01a2218e
sp_var = SingleVariable("u10")

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
ptbl_var = SingleVariable("utbl",path=srcdir())

# ╔═╡ 09204c11-38a7-47ea-b33a-b649dcac4d5b
pbl_var = SingleVariable("ubl",path=srcdir())

# ╔═╡ c4d1a045-ee7d-4235-9e84-1332b2616cff
blh_var = SingleVariable("blh")

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Timeseries Data for Stations
"

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
    ds = read_climatology(armsite,e5ds,ptbl_var)
    dtv = ds["valid_time"][:]
    ptbl = ds[ptbl_var.ncID][:]
    close(ds)
    ds = read_climatology(armsite,e5ds,sp_var)
    sp = ds[sp_var.ncID][:]
    close(ds)
    ds = read_climatology(armsite,e5ds,pbl_var)
    pbl = ds[pbl_var.ncID][:]
    close(ds)
    ds = read_climatology(armsite,e5ds,blh_var)
    blh = ds[blh_var.ncID][:]
    close(ds)
    md"Loading ERA5 Data"
end

# ╔═╡ 7cf5b637-8443-45b3-91db-781f01f75810
dtbeg = Date(2011,3)

# ╔═╡ 9b857d52-fab4-4b14-a58a-adcf9aee96b2
dtend = Date(2011,4)

# ╔═╡ c5bfe5af-2b16-4d5a-bf05-976b54915bea
begin
    ii = (dtv .>= dtbeg) .& (dtv .<= dtend)
    f1 = Figure()

    ax1_1 = Axis(
        f1[1,1],width=500,height=100,ylabel="Pressure / hPa",
        # yscale=log10
    )

    ax1_2 = Axis(
        f1[2,1],width=500,height=100,ylabel="Height / m",
    )

    lines!(ax1_1,dtv[ii],abs.(ptbl[ii]))
    lines!(ax1_1,dtv[ii],abs.(sp[ii]))
    lines!(ax1_1,dtv[ii],abs.(pbl[ii]))
    lines!(ax1_2,dtv[ii],blh[ii])
    # ylims!(ax1_1,1010,900)

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
# ╟─e282e373-f8a9-4a11-b1a7-e1d270b14090
# ╟─c3329018-a49d-4bc3-b6d3-5bc4ea401157
# ╠═bd0593a8-c898-4707-90cc-6ddb01a2218e
# ╠═aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╠═09204c11-38a7-47ea-b33a-b649dcac4d5b
# ╠═c4d1a045-ee7d-4235-9e84-1332b2616cff
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╠═a0757349-bf33-41bf-952e-723d874b9dcc
# ╠═7cf5b637-8443-45b3-91db-781f01f75810
# ╠═9b857d52-fab4-4b14-a58a-adcf9aee96b2
# ╠═c5bfe5af-2b16-4d5a-bf05-976b54915bea
