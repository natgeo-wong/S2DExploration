using DrWatson
@quickactivate "S2DExploration"
using ERA5Reanalysis
using DelimitedFiles
using Logging

include(srcdir("common.jl"))

e5ds = ERA5Hourly(start=Date(2001),stop=Date(2024,12,31),path=datadir())
evar = PressureVariable("t")

dtvec = e5ds.start : Day(1) : e5ds.stop; ndt = length(dtvec)
gID,lon,lat = readdlm("info.txt",',',comments=true,comment_char='#')
ipnt = closestnativelonlat(Point2(lon,lat))
vmat = zeros(37,24*ndt)

for idt in 1 : ndt
    @info "$(now()) - S2DExploration - Extraction point data $(evar.name) for ($(lon),$(lat)) from the DKRZ servers for $(dtvec[idt])"
    ibeg = (idt-1) * 24 + 1
    iend = idt * 24
    gds = dkrz(e5ds,evar,dtvec[idt])
    vmat[:,ibeg:iend] .= gds[evar.ID][ipnt,:]
    close(gds)
end

save_climatology(gID,e5ds,evar,vmat,ERA5Reanalysis.era5Pressures())