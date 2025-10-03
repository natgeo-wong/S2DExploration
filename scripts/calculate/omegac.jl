using DrWatson
@quickactivate "S2DExploration"

using ERA5Reanalysis

include(srcdir("calculate.jl"))

ID = "SGP"
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12,31),path=datadir())
calculate_ωc(e5ds,ID=ID)
calculate_ωc(e5ds,ID=ID,days=1)
calculate_ωc(e5ds,ID=ID,days=3)
calculate_ωc(e5ds,ID=ID,days=7)
calculate_ωc(e5ds,ID=ID,days=30)