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
# 04c. Wind Shear
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
ubl_var = SingleVariable("ubl",path=srcdir())

# ╔═╡ aab5e60b-6ffb-434c-a2d5-654cf9e5d082
vbl_var = SingleVariable("vbl",path=srcdir())

# ╔═╡ 09204c11-38a7-47ea-b33a-b649dcac4d5b
u_var = PressureVariable("u")

# ╔═╡ c4d1a045-ee7d-4235-9e84-1332b2616cff
v_var = PressureVariable("v")

# ╔═╡ 35916a1b-d7c5-4f5e-a68d-0bf2d21c9b87
ipre = ERA5Reanalysis.era5Pressures() .== 500

# ╔═╡ 59e95316-9d0c-462b-aedd-669ede68adcf
md"
### C. Timeseries Data for Stations
"

# ╔═╡ a0757349-bf33-41bf-952e-723d874b9dcc
begin
    ds = read_climatology("SGP",e5ds,ubl_var)
    dtv = ds["valid_time"][:]
    ubl = ds[ubl_var.ncID][:]
    close(ds)
    ds = read_climatology("SGP",e5ds,vbl_var)
    vbl = ds[vbl_var.ncID][:]
    close(ds)
    ds = read_climatology("SGP",e5ds,u_var)
    u500 = ds[u_var.ncID][ipre,:][:]
    close(ds)
    ds = read_climatology("SGP",e5ds,v_var)
    v500 = ds[v_var.ncID][ipre,:][:]
    close(ds)
    shear_SGP = sqrt.((ubl.-u500).^2 .+ (vbl.-v500).^2)
    ds = read_climatology("BNF",e5ds,ubl_var)
    dtv = ds["valid_time"][:]
    ubl = ds[ubl_var.ncID][:]
    close(ds)
    ds = read_climatology("BNF",e5ds,vbl_var)
    vbl = ds[vbl_var.ncID][:]
    close(ds)
    ds = read_climatology("BNF",e5ds,u_var)
    u500 = ds[u_var.ncID][ipre,:][:]
    close(ds)
    ds = read_climatology("BNF",e5ds,v_var)
    v500 = ds[v_var.ncID][ipre,:][:]
    close(ds)
    shear_BNF = sqrt.((ubl.-u500).^2 .+ (vbl.-v500).^2)
    md"Loading ERA5 Data"
end

# ╔═╡ c5bfe5af-2b16-4d5a-bf05-976b54915bea
begin
    ii = (month.(dtv).==4) .|| (month.(dtv).==5) .|| (month.(dtv).==6) .|| (month.(dtv).==7) .|| (month.(dtv).==8) .|| (month.(dtv).==9) .|| (month.(dtv).==10)
    # ii = (month.(dtv).==1) .|| (month.(dtv).==2) .|| (month.(dtv).==3) .|| (month.(dtv).==11) .|| (month.(dtv).==12)
    f1 = Figure()

    ax1_1 = Axis(
        f1[1,1],width=400,height=250,
        ylabel="Normalized Density",
        xlabel=L"Wind Shear Magnitude |$u_{500}-u_\text{bl}$| / m s$^{-1}$",xscale=log10
    )

    h1 = fit(Histogram,shear_SGP[ii],10. .^(-1:0.02:2),)
    h2 = fit(Histogram,shear_BNF[ii],10. .^(-1:0.02:2),)

    lines!(10. .^(-.99:0.02:2),h1.weights./sum(h1.weights)*150,label="SGP")
    lines!(10. .^(-.99:0.02:2),h2.weights./sum(h2.weights)*150,label="BNF")
    # ylims!(ax1_1,1010,900)

    axislegend(framevisible=false)

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
# ╟─bd0593a8-c898-4707-90cc-6ddb01a2218e
# ╟─aab5e60b-6ffb-434c-a2d5-654cf9e5d082
# ╟─09204c11-38a7-47ea-b33a-b649dcac4d5b
# ╟─c4d1a045-ee7d-4235-9e84-1332b2616cff
# ╟─35916a1b-d7c5-4f5e-a68d-0bf2d21c9b87
# ╟─59e95316-9d0c-462b-aedd-669ede68adcf
# ╟─a0757349-bf33-41bf-952e-723d874b9dcc
# ╟─c5bfe5af-2b16-4d5a-bf05-976b54915bea
