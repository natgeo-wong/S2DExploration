using DrWatson
@quickactivate "S2DExploration"

using RegionGrids, NASAMergedTb
using DelimitedFiles

geo = GeoRegion("SGP",path=srcdir())
tbd = TbDataset(start=Date(2001),stop=Date(2024,12,31),path=datadir())

sinfo = readdlm(datadir("ARMstations_$(geo.ID).csv"),',',skipstart=1)[:,[1,5,6]]
sID  = @views sinfo[:,1]

fol = joinpath(tbd.path,"timeseries-OLR")
if !isdir(fol); mkpath(joinpath(tbd.path,"timeseries-OLR")) end

for istn = 1 : nstn

    ds = NCDataset(joinpath(fol,"Tb-$(geo.ID)_$(sID[istn])-$(Dates.format(tbd.start,dateformat"yyyymmdd"))-$(Dates.format(tbd.stop,dateformat"yyyymmdd")).nc"))
    dt = ds["valid_time"].var[:]
    Tb = nomissing(ds["Tb"][:],NaN)
    close(ds)

    fnc = joinpath(tbd.path,"timeseries-OLR","OLR-$(geo.ID)_$(sID[istn])-$(Dates.format(tbd.start,dateformat"yyyymmdd"))-$(Dates.format(tbd.stop,dateformat"yyyymmdd")).nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now())",
    ))

    ds.dim["valid_time"] = ndt * 48

    nctime = defVar(ds,"valid_time",Float32,("valid_time",),attrib = Dict(
        "units"     => "hours since $(tbd.start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    ncvar = defVar(ds,"OLR",Float32,("valid_time",),attrib = Dict(
        "units"         => "K",
        "standard_name" => "brightness_temperature",
        "full_name"     => "Brightness Temperature",
    ))

    nctime.var[:] = dt
    ncvar[:] = 5.67e-8 * (Tb .* (1.228 .- 1.106e-3*Tb)).^4

    close(ds)

end