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

# ╔═╡ 16c7e940-f1b5-11ef-027b-ffed5bebb310
begin
    using Pkg; Pkg.activate()
    using DrWatson
end

# ╔═╡ 2dd9c606-8ab2-484e-8155-f2c908418fe0
begin
    @quickactivate "S2DExploration"
    using PlutoUI
    using Dates, DelimitedFiles
    using GeoRegions, ERA5Reanalysis, ETOPO
    using JSON3, Distances
    using CairoMakie, LaTeXStrings
    set_theme!(theme_latexfonts())
    md"Activating Project Environment for S2DExploration ..."
end

# ╔═╡ 2f6fac9c-a075-403b-a696-1fcb095f72d8
md"
# 01. Spatial Distribution of ARM Stations
"

# ╔═╡ a6877704-98ed-44fd-9ef7-55ec19202079
TableOfContents()

# ╔═╡ 13c98981-24fe-454b-8f46-28a243525cd9
md"
### A. Loading ARM Station and Plotting Data
"

# ╔═╡ 74a03763-63d0-4a16-9586-f55a5a3aed6f
geo_plot = GeoRegion([-101,-85,-85,-101,-101],[31,31,41,41,31])

# ╔═╡ 87a12c55-4dd4-4345-a534-d08b00b59bfc
geo_SGP = GeoRegion("plot_SGP",path=srcdir())

# ╔═╡ bfc8e1e2-f0d3-42c4-b7f2-3af00b747194
geo_BNF = GeoRegion("plot_BNF",path=srcdir())

# ╔═╡ 7038130a-cfcc-4e3b-a7e7-64157f0ab69b
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

# ╔═╡ 72ac56df-7c31-408a-85c4-7d7a4930676a
begin
    elon,elat = nativelonlat()
    md"Loading the native longitude/latitude points for the ERA5 Reanalysis model ..."
end

# ╔═╡ d3245994-1745-4f50-9512-14db7ad76527
begin
    states = JSON3.read(datadir("states.json"));
    md"Loading boundary data for states ..."
end

# ╔═╡ cd72582b-11cf-4697-887c-0b76fd88fe86
etpd = ETOPODataset(path=datadir())

# ╔═╡ 9e7eb0df-58b6-4bad-9e85-709647576360
lsd = getLandSea(etpd,geo_plot,save=false)

# ╔═╡ 838c2f2e-0a56-49db-9d50-992d44471cf0
lsd1 = getLandSea(etpd,GeoRegion([0,1,1,0,0].-87.75,[0,0,1,1,0].+33.9),save=false)

# ╔═╡ 5f910c13-5afe-49d7-b5e7-4bffba5b3ab4
lsd2 = getLandSea(etpd,GeoRegion([0,1,1,0,0]*2.5.-98.3,[0,0,1,1,0].+36.1),save=false)

# ╔═╡ 6bf4dda4-af40-4712-b290-9502202c70b3
lsd_SGP = getLandSea(etpd,geo_SGP)

# ╔═╡ dd2d3404-7cc5-40cf-a961-156598a6c38f
lsd_BNF = getLandSea(etpd,geo_BNF)

# ╔═╡ c8524322-c53d-485e-bd68-e154add31d3b
begin
    coast = readdlm(datadir("coast.cst"),comments=true)
    xc = coast[:,1]
    yc = coast[:,2]
    md"Loading Coastline Data ..."
end

# ╔═╡ 33d3451f-3299-4084-a986-74be3fe8ab2c
begin
    lon1,lat1 = coordinates(GeoRegion("AL",path=srcdir()));
    lon2,lat2 = coordinates(GeoRegion("LA",path=srcdir()));
    md"Loading state boundaries ..."
end

# ╔═╡ 9ec24adf-9adf-4f31-99b9-f7a13e1efa2d
md"
### B. Spatial Distribution of ARM Stations
"

# ╔═╡ a42c2eba-9abe-49a1-91db-b2c8eee8d804
begin
    valid1 = (sgp_info[:,3] .< Date(2001)) .& (sgp_info[:,4] .> Date(2005))
    valid2 = (sgp_info[:,3] .< Date(2021)) .& (sgp_info[:,4] .> Date(2025))
    md"Separating out the two distinct time periods of Southern Great Plains ARM Stations"
end

# ╔═╡ b33690e2-9200-4a00-b010-ca0fe45bfd49
begin
    fig = Figure()
    
    ax1 = Axis(
        fig[1,1],width=560,height=240,
        limits=(-100,-86,33,39),xticks=-100:5:-70,yticks=25:5:40,
        xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$"
    )
    heatmap!(ax1,lsd.lon,lsd.lat,lsd.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo,alpha=0.5)
    c = heatmap!(ax1,lsd1.lon,lsd1.lat,lsd1.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    heatmap!(ax1,lsd2.lon,lsd2.lat,lsd2.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax1,sgp_info[valid1,5],sgp_info[valid1,6],color=:blue)
    scatter!(ax1,sgp_info[valid2,5],sgp_info[valid2,6],color=:black)
    scatter!(ax1,sgp_info[1,5],sgp_info[1,6],color=:red,markersize=15)
    # lines!(ax1,xc,yc,color=:black)
    lines!(ax1,lon1,lat1,color=:black,linewidth=1)
    lines!(ax1,lon2,lat2,color=:black,linewidth=1)
    lines!(ax1,[-87.75,-92.5],[33.9,34.3],color=:black,linewidth=1,linestyle=(:dash,:dense))
    lines!(ax1,[-86.75,-88.13],[34.9,38.68],color=:black,linewidth=1,linestyle=(:dash,:dense))
    lines!(ax1,[-98.3,-99.57],[37.1,35.45],color=:black,linewidth=1,linestyle=(:dash,:dense))
    lines!(ax1,[-95.8,-93.38],[37.1,35.45],color=:black,linewidth=1,linestyle=(:dash,:dense))

    for istate in 1 : length(states.features)
        for ilines in 1 : length(states.features[istate].geometry.coordinates)
            coords = states.features[istate].geometry.coordinates[ilines]
            ncoord = length(coords); x = zeros(ncoord); y = zeros(ncoord)
            for icoord in 1 : ncoord
                x[icoord] = coords[icoord][1][1]
                y[icoord] = coords[icoord][2][1]
            end
            lines!(ax1,x,y,color=:black,linewidth=1)
        end
    end
    
    scatter!(ax1,bnf_info[:,5],bnf_info[:,6],color=:black)
    scatter!(ax1,bnf_info[1,5],bnf_info[1,6],color=:red,markersize=15)
    
    lines!(ax1,[0,1,1,0,0].-87.75,[0,0,1,1,0].+33.9,color=:black,linewidth=1)
    lines!(ax1,[0,1,1,0,0]*2.5.-98.3,[0,0,1,1,0].+36.1,color=:black,linewidth=1)

    ax_inset = Axis(
        fig[1, 1],width=175,height=175,halign=0.78,valign=0.8,
        limits=(-87.75,-86.75,33.9,34.9),
        xticksvisible=false,xticklabelsvisible=false,
        yticksvisible=false,yticklabelsvisible=false
    )

    heatmap!(ax_inset,lsd_BNF.lon,lsd_BNF.lat,lsd_BNF.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax_inset,bnf_info[1,5],bnf_info[1,6],color=:red,markersize=15)
    scatter!(ax_inset,bnf_info[2:end,5],bnf_info[2:end,6],color=:black)

    ax_inset = Axis(
        fig[1, 1],width=250,height=80,halign=0.05,valign=0.1,
        limits=(-98.3,-95.8,36.1,37.1),
        xticksvisible=false,xticklabelsvisible=false,
        yticksvisible=false,yticklabelsvisible=false
    )

    heatmap!(ax_inset,lsd_SGP.lon,lsd_SGP.lat,lsd_SGP.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax_inset,sgp_info[valid2,5],sgp_info[valid2,6],color=:black)
    scatter!(ax_inset,sgp_info[1,5],sgp_info[1,6],color=:red,markersize=15)

    Colorbar(fig[1,2], c, label="Topographic Height / km")
    resize_to_layout!(fig)
    fig
end

# ╔═╡ 3b93f80c-ddab-4c6c-ad79-feb894d9e6a0
CairoMakie.save(plotsdir("01-ARMstations.png"),fig)

# ╔═╡ d3b300e0-598b-4f8d-98f1-f69a8483e72b
md"
It can therefore be seen here that initially the ARM stations of the Southern Great Plains facility were much more widely spread out than those of Bankhead National Forest, which means that there is more opportunity there to study intra-grid variability for ERA5.

However, that doesn't mean that there aren't locations in the Southern Great Plains that can't do that, it just means that we cannot compare these locations to the Central Facility data.

Also, if you look at the ARM stations that are still available today in the Southern Great Plains, these ARM stations are distributed much more densely, and some can even be compared with the Central Facility, probably.
"

# ╔═╡ 3ead2b4a-bb47-4ecc-95cc-4959abe18823
md"
### C. Finding Neighbouring Stations (Intra ERA5-grid)?
"

# ╔═╡ 32ac38a0-036d-41b3-adef-30c7e018a696
begin
    sgp_info1 = sgp_info[valid1,:]; nsgp1 = size(sgp_info1,1)
    @bind iisgp1 Select([ii => "($(sgp_info1[ii,1])) $(sgp_info1[ii,2])" for ii in 1 : nsgp1])
end

# ╔═╡ 9b44c54b-ce1a-4ff8-9f69-70d6aceedd5f
begin
    sgp_info2 = sgp_info[valid2,:]; nsgp2 = size(sgp_info2,1)
    @bind iisgp2 Select([ii => "($(sgp_info2[ii,1])) $(sgp_info2[ii,2])" for ii in 1 : nsgp2])
end

# ╔═╡ 3dc8bd5f-18cf-4e24-90c4-a627b6035def
@bind iibnf Select([ii => "($(bnf_info[ii,1])) $(bnf_info[ii,2])" for ii in 1 : nbnf])

# ╔═╡ 727ead81-3dda-4d35-a1cb-817456c5bc05
begin
    sgp_dist1 = zeros(nsgp1,nsgp1)
    for istn in 1 : nsgp1, jstn in 1 : nsgp1
        sgp_dist1[jstn,istn] = haversine((sgp_info1[istn,5],sgp_info1[istn,6]),(sgp_info1[jstn,5],sgp_info1[jstn,6])) /1e3
    end
end

# ╔═╡ efbe37ab-63a6-44c8-80cf-dfa08d45e65e
begin
    sgp_dist2 = zeros(nsgp2,nsgp2)
    for istn in 1 : nsgp2, jstn in 1 : nsgp2
        sgp_dist2[jstn,istn] = haversine((sgp_info2[istn,5],sgp_info2[istn,6]),(sgp_info2[jstn,5],sgp_info2[jstn,6])) /1e3
    end
end

# ╔═╡ 59760189-d313-4987-b1c1-0ab77097a8f7
begin
    bnf_dist = zeros(nbnf,nbnf)
    for istn in 1 : nbnf, jstn in 1 : nbnf
        bnf_dist[jstn,istn] = haversine((bnf_info[istn,5],bnf_info[istn,6]),(bnf_info[jstn,5],bnf_info[jstn,6])) /1e3
    end
end

# ╔═╡ 8bb06a34-6999-4435-be1f-2fe2b2002d06
begin
    fig2 = Figure()
    
    ax2_1 = Axis(
        fig2[1,1],width=200,height=200,title="Southern Great Plains (Old)",
        limits=(-100.2,-94.8,33.8,39.2),xticks=-100:95,yticks=34:39
    )
    heatmap!(ax2_1,lsd_SGP.lon,lsd_SGP.lat,lsd_SGP.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax2_1,sgp_info1[:,5],sgp_info1[:,6],color=:black)
    scatter!(ax2_1,sgp_info1[sgp_dist1[iisgp1,:].<20,5],sgp_info1[sgp_dist1[iisgp1,:].<20,6],color=:blue)
    scatter!(ax2_1,sgp_info1[iisgp1,5],sgp_info1[iisgp1,6],color=:red,markersize=12)
    scatter!(ax2_1,elon.-360,elat,markersize=5,color=:grey)
    
    ax2_2 = Axis(
        fig2[1,2],width=200,height=200,title="Southern Great Plains (New)",
        limits=(-98.2,-95.8,35.3,37.7),xticks=-100:0.5:95,yticks=34:0.5:39
    )
    heatmap!(ax2_2,lsd_SGP.lon,lsd_SGP.lat,lsd_SGP.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax2_2,sgp_info2[:,5],sgp_info2[:,6],color=:black)
    scatter!(ax2_2,sgp_info2[sgp_dist2[iisgp2,:].<20,5],sgp_info2[sgp_dist2[iisgp2,:].<20,6],color=:blue)
    scatter!(ax2_2,sgp_info2[iisgp2,5],sgp_info2[iisgp2,6],color=:red,markersize=12)
    scatter!(ax2_2,elon.-360,elat,markersize=5,color=:grey)
    
    ax2_3 = Axis(
        fig2[1,3],width=200,height=200,title="Bankhead National Forest",
        limits=(-88.2,-85.8,33.3,35.7),xticks=-90:0.5:-85,yticks=32:0.5:37,
    )
    c2 = heatmap!(ax2_3,lsd_BNF.lon,lsd_BNF.lat,lsd_BNF.z ./1e3,colorrange=(-0.75,0.75),colormap=:topo)
    scatter!(ax2_3,bnf_info[:,5],bnf_info[:,6],color=:black)
    scatter!(ax2_3,bnf_info[bnf_dist[iibnf,:].<20,5],bnf_info[bnf_dist[iibnf,:].<20,6],color=:blue)
    scatter!(ax2_3,bnf_info[iibnf,5],bnf_info[iibnf,6],color=:red,markersize=12)
    scatter!(ax2_3,elon.-360,elat,markersize=5,color=:grey)
    
    Label(fig2[2,:],L"Longitude / $\degree$")
    Label(fig2[1,0],L"Latitude / $\degree$",rotation=pi/2)

    Colorbar(fig2[1,4], c2, label="Topographic Height / km")
    resize_to_layout!(fig2)
    fig2
end

# ╔═╡ f75b6f0f-fed5-46a2-ba79-fcda67637bb4
CairoMakie.save(plotsdir("01-nearbyARM.png"),fig2)

# ╔═╡ 179db92f-2ca8-4979-bc9f-aa613ae416fe
md"Save Specific Figure? $(@bind dosave PlutoUI.Slider(0:1))"

# ╔═╡ ae1b3b62-c2c8-4b0e-8b73-7051c4fd69be
if isone(dosave)
    CairoMakie.save(plotsdir("01-nearbyARM-$(sgp_info1[iisgp1,1])$(sgp_info2[iisgp2,1])$(bnf_info[iibnf,1]).png"),fig2)
end

# ╔═╡ e5a87bf6-09a5-420b-a971-51a278a11528
begin
    open(datadir("distances_sgpold.txt"),"w") do io
        writedlm(io,[sgp_info1[:,1] round.(sgp_dist1,digits=3)],',')
    end
    md"Saving distances for Old Southern Great Plains stations ..."
end

# ╔═╡ a125dca0-65d2-4e53-8eb5-71a65c4e4a90
begin
    open(datadir("distances_sgpnew.txt"),"w") do io
        writedlm(io,[sgp_info2[:,1] round.(sgp_dist2,digits=3)],',')
    end
    md"Saving distances for New Southern Great Plains stations ..."
end

# ╔═╡ ad899886-5ff2-4b8a-aa37-8544f185a675
begin
    open(datadir("distances_bnf.txt"),"w") do io
        writedlm(io,[bnf_info[:,1] round.(bnf_dist,digits=3)],',')
    end
    md"Saving distances for Bankhead National Forest stations ..."
end

# ╔═╡ Cell order:
# ╟─2f6fac9c-a075-403b-a696-1fcb095f72d8
# ╟─16c7e940-f1b5-11ef-027b-ffed5bebb310
# ╟─2dd9c606-8ab2-484e-8155-f2c908418fe0
# ╟─a6877704-98ed-44fd-9ef7-55ec19202079
# ╟─13c98981-24fe-454b-8f46-28a243525cd9
# ╟─74a03763-63d0-4a16-9586-f55a5a3aed6f
# ╟─87a12c55-4dd4-4345-a534-d08b00b59bfc
# ╟─bfc8e1e2-f0d3-42c4-b7f2-3af00b747194
# ╟─7038130a-cfcc-4e3b-a7e7-64157f0ab69b
# ╟─72ac56df-7c31-408a-85c4-7d7a4930676a
# ╟─d3245994-1745-4f50-9512-14db7ad76527
# ╟─cd72582b-11cf-4697-887c-0b76fd88fe86
# ╟─9e7eb0df-58b6-4bad-9e85-709647576360
# ╟─838c2f2e-0a56-49db-9d50-992d44471cf0
# ╟─5f910c13-5afe-49d7-b5e7-4bffba5b3ab4
# ╟─6bf4dda4-af40-4712-b290-9502202c70b3
# ╟─dd2d3404-7cc5-40cf-a961-156598a6c38f
# ╟─c8524322-c53d-485e-bd68-e154add31d3b
# ╟─33d3451f-3299-4084-a986-74be3fe8ab2c
# ╟─9ec24adf-9adf-4f31-99b9-f7a13e1efa2d
# ╟─a42c2eba-9abe-49a1-91db-b2c8eee8d804
# ╟─b33690e2-9200-4a00-b010-ca0fe45bfd49
# ╟─3b93f80c-ddab-4c6c-ad79-feb894d9e6a0
# ╟─d3b300e0-598b-4f8d-98f1-f69a8483e72b
# ╟─3ead2b4a-bb47-4ecc-95cc-4959abe18823
# ╟─32ac38a0-036d-41b3-adef-30c7e018a696
# ╟─9b44c54b-ce1a-4ff8-9f69-70d6aceedd5f
# ╟─3dc8bd5f-18cf-4e24-90c4-a627b6035def
# ╟─727ead81-3dda-4d35-a1cb-817456c5bc05
# ╟─efbe37ab-63a6-44c8-80cf-dfa08d45e65e
# ╟─59760189-d313-4987-b1c1-0ab77097a8f7
# ╟─8bb06a34-6999-4435-be1f-2fe2b2002d06
# ╟─f75b6f0f-fed5-46a2-ba79-fcda67637bb4
# ╟─179db92f-2ca8-4979-bc9f-aa613ae416fe
# ╟─ae1b3b62-c2c8-4b0e-8b73-7051c4fd69be
# ╟─e5a87bf6-09a5-420b-a971-51a278a11528
# ╟─a125dca0-65d2-4e53-8eb5-71a65c4e4a90
# ╟─ad899886-5ff2-4b8a-aa37-8544f185a675
