using DrWatson
@quickactivate "S2DExploration"
using ARMLive
using DelimitedFiles
using Logging
using Statistics

include(srcdir("common.jl"))

ads = ARMDataset(
    stream = "bnfdlprofwstats4newsM1.c1",
    start = Date(2024,9), stop = Date(2025,11,30), path = datadir()
)

mat = zeros(133,144,366)
count = zeros(133,144,366)
dtvec = ads.start : Day(1) : ads.stop
varID = "w"

for idt in dtvec
    try 
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
    catch
        @warn "reading file for $idt didn't work"
    end
end

mat ./= count

mat_monthly = zeros(133,144,12)

for imo = 1 : 12
    ndy = daysinmonth(Date(2004,imo))
    ii = dayofyear(Date(2004,imo,1)) : dayofyear(Date(2004,imo,ndy))
    matii = @views mat[:,:,ii]
    nlvl,nt,_ = size(matii)
    for it = 1 : nt, ilvl = 1 : nlvl
        iimat = @views matii[ilvl,it,:]
        mat_monthly[ilvl,it,imo] = mean(iimat[.!isnan.(iimat)])
    end
end

ds = read(ads,Date(2025),throw=false)
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