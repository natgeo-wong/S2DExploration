using Dates
using ERA5Reanalysis
using NCDatasets
using Printf

## DateString Aliasing
yrmo2dir(date::TimeType) = Dates.format(date,dateformat"yyyy/mm")
yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")
yr2str(date::TimeType)   = Dates.format(date,dateformat"yyyy")
ymd2str(date::TimeType)  = Dates.format(date,dateformat"yyyymmdd")

function read_climatology(
    armID :: String,
    e5ds  :: ERA5Hourly,
    evar  :: ERA5Variable;
    days  :: Int = 0
)

    if iszero(days)
        return NCDataset(joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * ".nc"
        ))
    else
        return NCDataset(joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        ))
    end

end

function save_climatology(
    armID :: String,
    e5ds  :: ERA5Hourly,
    evar  :: ERA5Variable,
    data  :: Vector{<:Real};
    days  :: Int = 0
)

    ndt = length(data)
    if iszero(days)
        fnc = joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * ".nc"
        )
    else
        fnc = joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        )
    end

    if isfile(fnc); rm(fnc) end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now())",
        "comments"    => "These NetCDF files were created in the same format that data is saved on the Climate Data Store",
    ))

    ds.dim["valid_time"] = ndt

    nctime = defVar(ds,"valid_time",Int64,("valid_time",),attrib = Dict(
        "units"     => "hours since $(e5ds.start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    ncvar = defVar(ds,evar.ID,Float64,("valid_time",),attrib = Dict(
        "long_name" => evar.long,
        "full_name" => evar.name,
        "units"     => evar.units,
    ))

    nctime[:] = collect(1:(ndt)) .- 1
    ncvar[:] = data

    close(ds)
    
end

function save_climatology(
    armID :: String,
    e5ds  :: ERA5Hourly,
    evar  :: ERA5Variable,
    data  :: Matrix{<:Real},
    lvls  :: Vector{Int};
    days  :: Int = 0
)

    nlvls,ndt = size(data)

    if iszero(days)
        fnc = joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * ".nc"
        )
    else
        fnc = joinpath(e5ds.path,"climatology",
            armID * "-" * evar.ID * "-" * 
            ymd2str(e5ds.start) * "-" * ymd2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        )
    end

    if isfile(fnc); rm(fnc) end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now())",
        "comments"    => "These NetCDF files were created in the same format that data is saved on the Climate Data Store",
    ))

    ds.dim["levels"]     = nlvls
    ds.dim["valid_time"] = ndt

    nctime = defVar(ds,"valid_time",Int64,("valid_time",),attrib = Dict(
        "units"     => "hours since $(e5ds.start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    nclvl = defVar(ds,"pressures",Int64,("levels",),attrib = Dict(
        "units"     => "hPa",
        "long_name" => "Pressure Levels",
    ))

    ncvar = defVar(ds,evar.ID,Float64,("levels","valid_time",),attrib = Dict(
        "long_name" => evar.long,
        "full_name" => evar.name,
        "units"     => evar.units,
    ))

    nctime[:]  = collect(1:ndt) .- 1
    nclvl[:]   = lvls
    ncvar[:,:] = data

    close(ds)
    
end