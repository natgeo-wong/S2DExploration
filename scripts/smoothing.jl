using DrWatson
@quickactivate "S2DExploration"

using RegionGrids, NASAMergedTb
using DelimitedFiles

include(srcdir("smoothing.jl"))

e5ds = ERA5Hourly(start=Date(2001),stop=Date(2024,12,31),path=datadir())
evar = SingleVariable("skt")

smoothing(e5ds,evar,ID="SGP",days=7)