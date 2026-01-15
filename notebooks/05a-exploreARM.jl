### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 5692ccd8-6056-11f0-07da-d1ae989cdac1
begin
    using Pkg; Pkg.activate()
    using DrWatson
end

# ╔═╡ bd5e50c0-5d29-410e-8389-beb8b636307d
begin
    @quickactivate "S2DExploration"
    using Dates
    using DelimitedFiles
    using ARMLive
    using Statistics
    
    using CairoMakie
    using LaTeXStrings
end

# ╔═╡ 618eccb3-457c-4530-8d33-4be497967800
md"
# 02. Explore Some Data
"

# ╔═╡ 43dd49a1-caef-4088-97a3-7ab2e4ae9569
ads = ARMDataset(
    stream="sgpinterpolatedsondeC1.c1",
    start=Date(2001),stop=Date(now()),path=datadir()
)

# ╔═╡ 91dae02d-9503-4a2b-8efc-a5db134a4237
dtii = Date(2016,4,5)

# ╔═╡ 4bedf8db-45f8-4332-a8df-bbbdf3426e31
begin
    ds = read(ads,dtii)
    t  = ds["time"][:]
    z  = ds["height"][:]
    t1 = nomissing(ds["temp"][:,:],NaN)
    sr = nomissing(ds["source_temp"][:,:],-1)
    # t2 = nomissing(ds["temp_sonde"][:,:],NaN)
    close(ds)
end

# ╔═╡ 76aca7c6-5b75-464a-94b0-1c3decd6d20b
ad2 = ARMDataset(
    stream="sgprlproftemp2news10mC1.c0",
    start=Date(2001),stop=Date(now()),path=datadir()
)

# ╔═╡ 26b596f2-1504-4c6b-83e7-f745e4c9504c
begin
    ds2 = read(ad2,dtii)
    dt2 = ds2["time"][:]
    z2  = ds2["height"][:]
    t2 = nomissing(ds2["temp_sonde"][:,:],NaN)
    close(ds2)
end

# ╔═╡ 5e3b22c2-e8f0-4735-b35f-ce9a9aa00bce
f1 = Figure();

# ╔═╡ ecd1ec06-6957-4c75-9a96-8a82bcd50598
begin
    empty!(f1)
    ax1 = Axis(f1[1,1],height=250,width=300,limits=(0,1440,0,18))
    
    contourf!(ax1,1:length(t),z,t1' .+273.15,levels=200:5:300)
    
    ax2 = Axis(f1[1,2],height=250,width=300,limits=(0,1440,0,18))
    
    contourf!(ax2,1:length(t),z,sr')
    
    ax3 = Axis(f1[1,3],height=250,width=300,limits=(0,1440,0,18))
    
    contourf!(ax3,((1:length(dt2)).-0.5)*10,z2,t2',levels=200:5:300)
    
    resize_to_layout!(f1)
    f1
end

# ╔═╡ 3690323e-8928-45d5-8334-0d20614193c0
t1

# ╔═╡ 1f072f07-7ace-43ac-a999-de4cce3b68a7
length(z2.<18)

# ╔═╡ Cell order:
# ╟─618eccb3-457c-4530-8d33-4be497967800
# ╟─5692ccd8-6056-11f0-07da-d1ae989cdac1
# ╟─bd5e50c0-5d29-410e-8389-beb8b636307d
# ╠═43dd49a1-caef-4088-97a3-7ab2e4ae9569
# ╠═91dae02d-9503-4a2b-8efc-a5db134a4237
# ╟─4bedf8db-45f8-4332-a8df-bbbdf3426e31
# ╠═76aca7c6-5b75-464a-94b0-1c3decd6d20b
# ╟─26b596f2-1504-4c6b-83e7-f745e4c9504c
# ╠═5e3b22c2-e8f0-4735-b35f-ce9a9aa00bce
# ╠═ecd1ec06-6957-4c75-9a96-8a82bcd50598
# ╠═3690323e-8928-45d5-8334-0d20614193c0
# ╠═1f072f07-7ace-43ac-a999-de4cce3b68a7
