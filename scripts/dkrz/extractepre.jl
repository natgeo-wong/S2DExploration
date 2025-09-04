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

fnc = gID * "-" * evar.ID * "-" * ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * ".nc"
fnc = joinpath(e5ds.path,fnc)

if isfile(fnc); rm(fnc) end
ds = NCDataset(fnc,"c",attrib = Dict(
    "Conventions" => "CF-1.6",
    "history"     => "Created on $(Dates.now()) with ERA5Reanalysis.jl",
    "comments"    => "ERA5Reanalysis.jl creates NetCDF files in the same format that data is saved on the Climate Data Store",
    "doi"         => e5ds.pldoi
))

ds.dim["levels"] = 37
ds.dim["valid_time"] = ndt * 24

nctime = defVar(ds,"valid_time",Int64,("valid_time",),attrib = Dict(
    "units"     => "hours since $(e5ds.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

nclvl = defVar(ds,"pressures",Int64,("levels",),attrib = Dict(
    "units"     => "hPa",
    "long_name" => "Pressure Levels",
))

ncvar = defVar(ds,evar.ID,Float64,("levels","valid_time"),attrib = Dict(
    "long_name"     => evar.long,
    "full_name"     => evar.name,
    "units"         => evar.units,
))

nctime[:] = collect(1:(24*ndt)) .- 1
nclvl[:] = ERA5Reanalysis.era5Pressures()
ncvar[:,:] = vmat

close(ds)