using DrWatson
@quickactivate "S2DExploration"

using RegionGrids, NASAMergedTb
using DelimitedFiles

include(srcdir("calculate.jl"))

geo = GeoRegion("SGP",path=srcdir())
tbd = TbDataset(start=Date(2001),stop=Date(2024,12,31),path=datadir())

sinfo = readdlm(datadir("ARMstations_$(geo.ID).csv"),',',skipstart=1)[1,[1,5,6]]
ipnt = closestnativelonlat(Point2(sinfo[2],sinfo[3]))

dtvec = tbd.start : Day(1) : tbd.stop; ndt = length(dtvec)

lon,lat = NASAMergedTb.btdlonlat()
ggrd = RegionGrid(geo,lon,lat)
nlon = length(ggrd.lon)
nlat = length(ggrd.lat)
blongrid = zeros(nlon,nlat); blongrid .= ggrd.lon
blatgrid = zeros(nlon,nlat); blatgrid .= reshape(ggrd.lat,1,:)

iarg = closestnativelonlat(Point2.(blongrid[:],blatgrid[:])) .== ipnt
tbdata = zeros(Float32,sum(iarg),48,ndt) * NaN

tmpd = zeros(Float32,nlon,nlat,48)
for idt = 1 : ndt

    ids = read(tbd,geo,dtvec[idt])
    NCDatasets.load!(ids["Tb"].var,tmpd,:,:,:)
    tbdata[:,:,idt] = reshape(tmpd,nlon*nlat,:)[iarg,:]
    close(ids)

end

fol = joinpath(tbd.path,"timeseries-Tb_compare")
if !isdir(fol); mkpath(joinpath(tbd.path,"timeseries-Tb_compare")) end

fnc = joinpath(fol,"Tb_compare-$(geo.ID)_$(sinfo[1])-$(Dates.format(tbd.start,dateformat"yyyymmdd"))-$(Dates.format(tbd.stop,dateformat"yyyymmdd")).nc")
if isfile(fnc); rm(fnc,force=true) end

ds = NCDataset(fnc,"c",attrib = Dict(
    "Conventions" => "CF-1.6",
    "history"     => "Created on $(Dates.now())",
))

ds.dim["coordinates"] = sum(iarg)
ds.dim["valid_time"] = ndt * 48

nctime = defVar(ds,"valid_time",Float32,("valid_time",),attrib = Dict(
    "units"     => "hours since $(tbd.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

nclon = defVar(ds,"longitude",Float32,("coordinates",),attrib = Dict(
    "units"     => "degrees_east",
    "long_name" => "longitude",
))

nclat = defVar(ds,"latitude",Float32,("coordinates",),attrib = Dict(
    "units"     => "degrees_north",
    "long_name" => "latitude",
))

ncvar = defVar(ds,"Tb",Float32,("coordinates","valid_time"),attrib = Dict(
    "units"         => "K",
    "standard_name" => "brightness_temperature",
    "full_name"     => "Brightness Temperature",
))

nctime[:] = (collect(1:(48*ndt)) .- 0.5) ./ 2
nclon[:] = blongrid[iarg]
nclat[:] = blatgrid[iarg]
ncvar[:] = reshape(tbdata,sum(iarg),:)

close(ds)