using DrWatson
@quickactivate "S2DExploration"
using ARMLive
using DelimitedFiles
using Logging
using Statistics

include(srcdir("common.jl"))

ads = ARMDataset(
    stream = "sgpdlprofwstats4newsC1.c1",
    start = Date(2000), stop = Date(2024,12,31), path = datadir()
)

mat = zeros(133,144,366)
count = zeros(133,144,366)
dtvec = ads.start : Day(1) : ads.stop
varID = "w"

for idt in dtvec
    ids = read(ads,idt,throw=false)
    if !isnothing(ids)
        if !isleapyear(idt)
            if month(idt) > 2
                ii = dayofyear(idt) + 1
            else
                ii = dayofyear(idt)
            end
        else
            ii = dayofyear(idt)
        end
        mat[:,:,ii] .+= nomissing(ids[varID][:,:],0)
        count[:,:,ii] += .!isnan.(nomissing(ids[varID][:,:],NaN))
        close(ids)
    end
end

mat ./= count

mat_monthly = zeros(133,1440,12)

for imo = 1 : 12
    ndy = daysinmonth(Date(2004,imo))
    ii = dayofyear(Date(2004,imo,1)) : dayofyear(Date(2004,imo,ndy))
    matii = @views mat[:,:,ii]
    if imo != 2
        mat_monthly[:,:,imo] = mean(matii,dims=3)
    else
        mat_monthly[:,:,imo] = sum(matii,dims=3) ./ 28.24 # Since 2000 is not a leap year
    end
end

ds = read(ads,Date(2024),throw=false)
z  = ds["height"][:]
varDict = Dict(ds[varID].attrib)
close(ds)

mat_monthly = reshape(mat_monthly,332,6,24,12)
mat_monthly = dropdims(mean(mat_monthly,dims=2),dims=2)

fnc = joinpath(
    ads.path,
    "monthhour-$(varID)-$(ymd2str(ads.start))_$(ymd2str(ads.stop)).nc"
)
if isfile(fnc); rm(fnc,force=true) end

ds = NCDataset(fnc,"c",attrib = Dict(
    "Conventions" => "CF-1.6",
    "history"     => "Created on $(Dates.now())"
))

ds.dim["month"] = 12
ds.dim["valid_time"] = 24
ds.dim["levels"] = 133

ncz = defVar(ds,"height",Float64,("levels",),attrib = Dict(
    "units"     => "km",
    "long_name" => "height",
))

ncvar = defVar(ds,varID,Float64,("levels","valid_time","month",),attrib = varDict)

ncz[:] = z
ncvar[:,:,:] = mat_monthly

close(ds)