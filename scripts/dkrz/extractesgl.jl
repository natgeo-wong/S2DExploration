using DrWatson
@quickactivate "S2DExploration"
using ERA5Reanalysis
using DelimitedFiles
using Logging

include(srcdir("common.jl"))

e5ds = ERA5Hourly(start=Date(2001),stop=Date(2024,12,31),path=datadir())
esgllist = [
    SingleVariable("tp"), SingleVariable("tcc"), SingleVariable("tcw"),
    SingleVariable("hcc"), SingleVariable("mcc"), SingleVariable("lcc"), 
    SingleVariable("t2m"), SingleVariable("d2m"), SingleVariable("skt"), 
    SingleVariable("blh"), SingleVariable("u10"), SingleVariable("v10"),
    SingleVariable("tciw"), SingleVariable("tclw"), SingleVariable("tcwv"),
    SingleVariable("sp"), SingleVariable("sshf"), SingleVariable("slhf"),
]

dtvec = e5ds.start : Day(1) : e5ds.stop; ndt = length(dtvec)
gID,lon,lat = readdlm("info.txt",',',comments=true,comment_char='#')
ipnt = closestnativelonlat(Point2(lon,lat))
vmat = zeros(24*ndt)

for evar in esgllist

    for idt in 1 : ndt
        @info "$(now()) - S2DExploration - Extraction point data $(evar.name) for ($(lon),$(lat)) from the DKRZ servers for $(dtvec[idt])"
        ibeg = (idt-1) * 24 + 1
        iend = idt * 24
        gds = dkrz(e5ds,evar,dtvec[idt])
        vmat[ibeg:iend] .= nomissing(gds[evar.ID][ipnt,:],NaN)
        close(gds)
    end

    save_climatology(gID,e5ds,evar,vmat)

end