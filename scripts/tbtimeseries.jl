using DrWatson
@quickactivate "S2DExploration"

using RegionGrids, NASAMergedTb
using DelimitedFiles

geo = GeoRegion("KS_OK",path=srcdir())
tbd = TbDataset(start=Date(2001),stop=Date(2024,12,31),path=datadir())

dtvec = tbd.start : Day(1) : tbd.stop; ndt = length(dtvec)
tbdata = zeros(Float32,48,ndt)

gID,slon,slat = readdlm("info.txt",',',comments=true,comment_char='#')

lon,lat = NASAMergedTb.btdlonlat()
ggrd = RegionGrid(geo,lon,lat)
ilon = argmin(abs.(slon.-ggrd.lon))
ilat = argmin(abs.(slat.-ggrd.lat))

for idt = 1 : ndt

    ids = read(tbd,geo,dtvec[idt])
    tbdata[:,idt] = nomissing(ids["Tb"][ilon,ilat,:],NaN)
    close(ids)

end

fnc = joinpath(tbd.path,"$gID-Tb-$(Dates.format(tbd.start,dateformat"yyyymmdd"))-$(Dates.format(tbd.stop,dateformat"yyyymmdd")).nc")
if isfile(fnc); rm(fnc,force=true) end

ds = NCDataset(fnc,"c",attrib = Dict(
    "Conventions" => "CF-1.6",
    "history"     => "Created on $(Dates.now())",
))

ds.dim["valid_time"] = ndt * 48

nctime = defVar(ds,"valid_time",Int64,("valid_time",),attrib = Dict(
    "units"     => "hours since $(tbd.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncvar = defVar(ds,"Tb",Float32,("valid_time",),attrib = Dict(
    "units"         => "K",
    "standard_name" => "brightness_temperature",
    "full_name"     => "Brightness Temperature",
))

nctime[:] = (collect(1:(48*ndt)) .- 1) ./ 2
ncvar[:] = tbdata[:]

close(ds)