using DrWatson
@quickactivate "S2DExploration"

using RegionGrids, NASAMergedTb
using DelimitedFiles

geo = GeoRegion("SGP",path=srcdir())
tbd = TbDataset(start=Date(2001),stop=Date(2024,12,31),path=datadir())

sinfo = readdlm(datadir("ARMstations_$(geo.ID).csv"),',',skipstart=1)[:,[1,5,6]]
sID  = @views sinfo[:,1]
slon = @views sinfo[:,2]
slat = @views sinfo[:,3]
nstn = length(sID)

dtvec = tbd.start : Day(1) : tbd.stop; ndt = length(dtvec)
tbdata = zeros(Float32,48,ndt,nstn) * NaN

lon,lat = NASAMergedTb.btdlonlat()
ggrd = RegionGrid(geo,lon,lat)
nlon = length(ggrd.lon)
nlat = length(ggrd.lat)

ilon = zeros(Int,nstn)
ilat = zeros(Int,nstn)
tmpd = zeros(Float32,nlon,nlat,48)

for istn = 1 : nstn
    if .!isnan(slon[istn]) && .!isnan(slat[istn])
        ilon[istn] = argmin(abs.(slon[istn].-ggrd.lon))
        ilat[istn] = argmin(abs.(slat[istn].-ggrd.lat))
    end
end

for idt = 1 : ndt

    ids = read(tbd,geo,dtvec[idt])
    NCDatasets.load!(ids["Tb"].var,tmpd,:,:,:)

    for istn = 1 : nstn
        if !iszero(ilon[istn])
            tbdata[:,idt,istn] = tmpd[ilon[istn],ilat[istn],:]
        end
    end
    
    close(ids)

end

for istn = 1 : nstn

    fnc = joinpath(tbd.path,"$(geo.ID)_$(sID[istn])-Tb-$(Dates.format(tbd.start,dateformat"yyyymmdd"))-$(Dates.format(tbd.stop,dateformat"yyyymmdd")).nc")
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

    ncvar = defVar(ds,"Tb",Float32,("valid_time",),attrib = Dict(
        "units"         => "K",
        "standard_name" => "brightness_temperature",
        "full_name"     => "Brightness Temperature",
    ))

    nctime[:] = (collect(1:(48*ndt)) .- 1) ./ 2
    ncvar[:] = tbdata[:,:,istn][:]

    close(ds)

end