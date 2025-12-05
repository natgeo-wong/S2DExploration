using DrWatson
@quickactivate "S2DExploration"
using ARMLive
using DelimitedFiles
using Logging
using Statistics

include(srcdir("common.jl"))

ads = ARMDataset(
    stream = "sgparmbeatmC1.c1",
    start = Date(2000), stop = Date(2023,12,31), path = datadir()
)

dtvec = ads.start : Year(1) : ads.stop
varID = "temp"
mat = []

for idt in dtvec
    ids = read(ads,idt,throw=false)
    mat = cat(mat,nomissing(ids[varID][:,:],NaN),dims=2)
end

ds = read(ads,Date(2000),throw=false)
varDict = Dict(ds[varID].attrib)
close(ds)

fnc = joinpath(
    ads.path,
    "timeseries-$(varID)-$(ymd2str(ads.start))_$(ymd2str(ads.stop)).nc"
)
if isfile(fnc); rm(fnc,force=true) end

ds = NCDataset(fnc,"c",attrib = Dict(
    "Conventions" => "CF-1.6",
    "history"     => "Created on $(Dates.now())"
))

ds.dim["valid_time"] = size(mat,2)
ds.dim["levels"] = 37

ncvar = defVar(ds,varID,Float64,("levels","valid_time"),attrib = varDict)
ncvar[:,:] = mat

close(ds)