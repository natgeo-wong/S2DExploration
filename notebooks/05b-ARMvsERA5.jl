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
    using Dates, DelimitedFiles, StatsBase
    using GeoRegions, ERA5Reanalysis, ARMLive
    using Dierckx
    using DataInterpolations
    using CairoMakie, LaTeXStrings
    set_theme!(theme_latexfonts())

    include(srcdir("smoothing.jl"))
    md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 05b. ARM vs ERA5 Comparison
"

# ╔═╡ 62a477f5-24a0-46d8-b861-cd8d5363fdc0
TableOfContents()

# ╔═╡ 1e8da1c7-591e-4a16-8d42-ca0bc8a26a59
md"
### A. Loading Station Information
"

# ╔═╡ 61f81912-baff-43cd-8f67-8993ad9a7429
@bind armsite Select([
    "SGP" => "(SGP) Southern Great Plains",
    "BNF" => "(BNF) Bankhead National Forest",
])

# ╔═╡ 1df29268-e315-48bd-ba69-9bd40809a6d0
getarmcoordinates() = if armsite == "SGP"
    return -97.487643, 36.607322;
else
    return -87.338177, 34.342481;
end

# ╔═╡ 587abd2c-301c-4fca-acf3-bcb3b10c8394
getarmdatastream() = if armsite == "SGP"
    return "sgpinterpolatedsondeC1.c1";
else
    return "bnfinterpolatedsondeM1.c1";
end

# ╔═╡ 848bee35-e9b0-4eb9-a367-6be0a7b579d6
getarmdldatastream() = if armsite == "SGP"
    return "sgpdlprofwstats4newsC1.c1";
else
    return "bnfdlprofwstats4newsM1.c1";
end

# ╔═╡ beef2002-39ee-462f-85d7-dda304648cac
getarmvaliddates() = if armsite == "SGP"
    return "20000101_20241231";
else
    return "20240901_20251130";
end

# ╔═╡ c5cab2fa-847f-4187-9d21-9fa864a8c3ff
md"
### B. Defining ERA5 Datasets and Variables
"

# ╔═╡ d6b40589-4329-479a-a4b9-3355b1731d12
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12),path=datadir())

# ╔═╡ 94c56035-1215-40a6-872e-42e15dcdf262
md"
### C. Defining Interpolation Functions to Local Time
"

# ╔═╡ 1330b240-d031-48d9-b3f4-7d6f34d2508d
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ 399dac03-a75b-4ab1-9801-5a91fa8f5420
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

# ╔═╡ 1bdcde35-336e-474c-a936-e761ded97893
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

# ╔═╡ d9f8c27a-2de4-49c8-ad7c-f28afeac2017
md"
### D. Loading Timeseries Data for Stations
"

# ╔═╡ 67171dd2-cac9-440b-9b2e-45533dab368a
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

# ╔═╡ dcff524d-19c3-4643-bba1-f8815fd05d52
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

# ╔═╡ a9aab853-930d-4565-8661-98e3e2e0593f
begin
    t_s,sp = load_climatology(armsite,e5ds,SingleVariable("sp"))
    μsp = mean(sp) ./ 100
    md"Loading single level variable climatologies ..."
end

# ╔═╡ f09d038e-803f-4270-89f6-8ddb07b2d63f
begin
    t_p,lvl,w = load_climatology(armsite,e5ds,PressureVariable("w"));
    _,_,t = load_climatology(armsite,e5ds,PressureVariable("t"));
    _,_,q = load_climatology(armsite,e5ds,PressureVariable("q"));
    md"Loading pressure level variable climatologies ..."
end

# ╔═╡ f7ef59db-5722-4acf-9574-ed25e7075e64
md"
### E. Loading ARM Observational Dataset for $(armsite)
"

# ╔═╡ 43dd49a1-caef-4088-97a3-7ab2e4ae9569
begin
    ads = ARMDataset(
        stream=getarmdatastream(),
        start=Date(2001),stop=Date(now()),path=datadir()
    )
    md"Loading the ARM Interpolated Sonde Dataset"
end

# ╔═╡ 4bedf8db-45f8-4332-a8df-bbbdf3426e31
begin
    ds  = NCDataset(joinpath(ads.path,"monthhour-bar_pres-$(getarmvaliddates()).nc"))
    z   = ds["height"][:]
    pre = ds["bar_pres"][:,:,:] * 10; μp = dropdims(mean(pre,dims=3),dims=3)
    p   = dropdims(mean(μp,dims=2),dims=2)
    close(ds)
    ds  = NCDataset(joinpath(ads.path,"monthhour-temp-$(getarmvaliddates()).nc"))
    ta  = ds["temp"][:,:,:]; μta = dropdims(mean(ta,dims=3),dims=3)
    close(ds)
    ds  = NCDataset(joinpath(ads.path,"monthhour-vap_pres-$(getarmvaliddates()).nc"))
    vp  = ds["vap_pres"][:,:,:] * 10; μvp = dropdims(mean(vp,dims=3),dims=3)
    close(ds)
    md"Loading compiled climatology ..."
end

# ╔═╡ 124b956d-4289-4035-83fe-2fa86ca1a20e
begin
    ads2 = ARMDataset(
        stream=getarmdldatastream(),
        start=Date(2001),stop=Date(now()),path=datadir()
    )
    md"Loading the ARM Doppler Lidar Vertical Velocity Dataset"
end

# ╔═╡ 6d51ba77-c8a0-43a9-91bd-f6f164fc15f2
begin
    ds2 = NCDataset(joinpath(ads2.path,"monthhour-w-$(getarmvaliddates()).nc"))
    z2  = ds2["height"][:]
    ω   = ds2["w"][:,:,:]; μω = dropdims(mean(ω,dims=3),dims=3)
    close(ds2)
    md"Loading compiled climatology ..."
end

# ╔═╡ ebdac4f1-4050-4cbb-9a95-a0da38fbbee6
begin
    μq = 0.622 .* μvp ./ (μp .- 0.378 .* μvp);
    md"Calculating the Specific Humidity"
end

# ╔═╡ 60d3fa1d-fc94-400c-b15b-15080d8c77a1
md"
### F. Converting $w$ to $\omega$ and z-coordinates to p-coordinates 
"

# ╔═╡ 3d33c631-f2e7-4597-b046-e7c630c76617
spl = Spline1D(z[.!isnan.(p)],p[.!isnan.(p)],k=2,bc="extrapolate"); p2 = spl.(z2/1000);

# ╔═╡ 6c0c991a-5ec6-47f8-937b-50db0121b7e4
μω2 = μω ./ (-9.81 .* p2 / 1000);

# ╔═╡ 8a03bc44-ffa5-474b-8462-f764d088e036
md"
### G. Figures!
"

# ╔═╡ 5e3b22c2-e8f0-4735-b35f-ce9a9aa00bce
fig = Figure();

# ╔═╡ ecd1ec06-6957-4c75-9a96-8a82bcd50598
begin
    empty!(fig)
    axs = Matrix{Axis}(undef,3,3)

    for irow = 1 : 3, icol = 2 : 3
        if isone(irow)
            axs[irow,icol] = Axis(
                fig[irow,icol],width=200,height=120,
                xticklabelsvisible=false,
                yticklabelsvisible=isone(icol),
                # yscale=log10,yticks=[10,100,1000]
            )
            ylims!(axs[irow,icol],1000,10)
        else
            axs[irow,icol] = Axis(
                fig[irow,icol],width=200,height=120,
                xticklabelsvisible=irow==3,
                yticklabelsvisible=isone(icol),
            )
            ylims!(axs[irow,icol],1000,10)
        end
        xlims!(axs[irow,icol],0,24)
    end

    for irow = 1 : 3
        axs[irow,1] = Axis(fig[irow,1],width=40,height=120)
        ylims!(axs[irow,1],1000,10)
    end

    lines!(axs[1,1],dropdims(mean(w,dims=2),dims=2),lvl,linewidth=3)
    lines!(axs[1,1],dropdims(mean(μω2[z2.<1500,:],dims=2),dims=2),p2[z2.<1500])
    xlims!(axs[1,1],-0.11,0.11)
    axs[1,1].xticks = [-0.1,0,0.1]
    lines!(axs[2,1],dropdims(mean(t,dims=2),dims=2),lvl,linewidth=3)
    lines!(axs[2,1],dropdims(mean(μta.+273.15,dims=2),dims=2),p)
    xlims!(axs[2,1],190,310)
    axs[2,1].xticks = [200,300]
    lines!(axs[3,1],dropdims(mean(q,dims=2),dims=2),lvl,linewidth=3)
    lines!(axs[3,1],dropdims(mean(μq,dims=2),dims=2),p)
    axs[3,1].xscale=log10

    cbr_1 = heatmap!(axs[1,2],t_p,lvl,w',colorrange=(-1,1).*1e-1,colormap=:RdBu)
    cbr_2 = heatmap!(axs[2,2],t_p,lvl,(t.-mean(t,dims=2))',colorrange=(-2,2),colormap=:balance)
    cbr_3 = heatmap!(axs[3,2],t_p,lvl,(q./mean(q,dims=2) .-1)'*100,colorrange=(-5,5).*2,colormap=:BrBg)

    heatmap!(axs[1,3],(1:24).-6.5 ,p2[z2.<1500],μω2[z2.<1500,:]',colorrange=(-1,1).*1e-1,colormap=:RdBu)
    heatmap!(axs[1,3],(1:24).+17.5,p2[z2.<1500],μω2[z2.<1500,:]',colorrange=(-1,1).*1e-1,colormap=:RdBu)
    
    heatmap!(axs[2,3],(1:24).-6.5,p,(μta.-mean(μta,dims=2))',colorrange=(-2,2),colormap=:balance)
    heatmap!(axs[2,3],(1:24).+17.5,p,(μta.-mean(μta,dims=2))',colorrange=(-2,2),colormap=:balance)

    heatmap!(axs[3,3],(1:24).-6.5,p,(μq./mean(μq,dims=2) .-1)'*100,colorrange=(-5,5).*2,colormap=:BrBg)
    heatmap!(axs[3,3],(1:24).+17.5,p,(μq./mean(μq,dims=2) .-1)'*100,colorrange=(-5,5).*2,colormap=:BrBg)

    Label(fig[1:3,0],"Pressure / hPa",rotation=pi/2)
    Label(fig[4,2:3],"Hour of Day")
    Label(fig[1,2],halign=:left,valign=:top,padding=(7,0,0,5),L"ERA5 $\omega$")
    Label(fig[2,2],halign=:left,valign=:top,padding=(7,0,0,5),L"ERA5 $T'$")
    Label(fig[3,2],halign=:left,valign=:top,padding=(7,0,0,5),L"ERA5 $\frac{q'}{\mu(q)} - 1$")
    Label(fig[1,3],halign=:left,valign=:top,padding=(7,0,0,5),L"Doppler Lidar $\omega$")
    Label(fig[2,3],halign=:left,valign=:top,padding=(7,0,0,5),L"Sonde $T'$")
    Label(fig[3,3],halign=:left,valign=:top,padding=(7,0,0,5),L"Sonde $\frac{q'}{\mu(q)} - 1$")

    Colorbar(fig[1,4],cbr_1,label=L"Pa s$^{-1}$")
    Colorbar(fig[2,4],cbr_2,label="K")
    Colorbar(fig[3,4],cbr_3,label=L"%")
    
    Label(fig[0,1:3],"ERA5 vs ARM Observations ($(armsite))",font=:bold)
    resize_to_layout!(fig)
    fig
end

# ╔═╡ d7de8970-0d83-4036-9c86-9b81c3e5b66c
CairoMakie.save(plotsdir("05b-ERA5vsARM-$armsite.png"),fig)

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╠═bd5e50c0-5d29-410e-8389-beb8b636307d
# ╟─62a477f5-24a0-46d8-b861-cd8d5363fdc0
# ╟─1e8da1c7-591e-4a16-8d42-ca0bc8a26a59
# ╟─61f81912-baff-43cd-8f67-8993ad9a7429
# ╟─1df29268-e315-48bd-ba69-9bd40809a6d0
# ╠═587abd2c-301c-4fca-acf3-bcb3b10c8394
# ╟─848bee35-e9b0-4eb9-a367-6be0a7b579d6
# ╠═beef2002-39ee-462f-85d7-dda304648cac
# ╟─c5cab2fa-847f-4187-9d21-9fa864a8c3ff
# ╟─d6b40589-4329-479a-a4b9-3355b1731d12
# ╟─94c56035-1215-40a6-872e-42e15dcdf262
# ╟─1330b240-d031-48d9-b3f4-7d6f34d2508d
# ╟─399dac03-a75b-4ab1-9801-5a91fa8f5420
# ╟─1bdcde35-336e-474c-a936-e761ded97893
# ╟─d9f8c27a-2de4-49c8-ad7c-f28afeac2017
# ╟─67171dd2-cac9-440b-9b2e-45533dab368a
# ╠═dcff524d-19c3-4643-bba1-f8815fd05d52
# ╠═a9aab853-930d-4565-8661-98e3e2e0593f
# ╠═f09d038e-803f-4270-89f6-8ddb07b2d63f
# ╟─f7ef59db-5722-4acf-9574-ed25e7075e64
# ╠═43dd49a1-caef-4088-97a3-7ab2e4ae9569
# ╠═4bedf8db-45f8-4332-a8df-bbbdf3426e31
# ╠═124b956d-4289-4035-83fe-2fa86ca1a20e
# ╠═6d51ba77-c8a0-43a9-91bd-f6f164fc15f2
# ╟─ebdac4f1-4050-4cbb-9a95-a0da38fbbee6
# ╟─60d3fa1d-fc94-400c-b15b-15080d8c77a1
# ╠═3d33c631-f2e7-4597-b046-e7c630c76617
# ╠═6c0c991a-5ec6-47f8-937b-50db0121b7e4
# ╟─8a03bc44-ffa5-474b-8462-f764d088e036
# ╟─5e3b22c2-e8f0-4735-b35f-ce9a9aa00bce
# ╟─ecd1ec06-6957-4c75-9a96-8a82bcd50598
# ╟─d7de8970-0d83-4036-9c86-9b81c3e5b66c
