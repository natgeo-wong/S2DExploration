using DrWatson
@quickactivate "S2DExploration"

using ERA5Reanalysis

include(srcdir("calculate.jl"))

ID = "SGP"
e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12,31),path=datadir())

calculate_utbl(e5ds,ID=ID)
calculate_utbl(e5ds,ID=ID,days=1)
calculate_utbl(e5ds,ID=ID,days=3)
calculate_utbl(e5ds,ID=ID,days=7)
calculate_utbl(e5ds,ID=ID,days=30)

calculate_ubl(e5ds,ID=ID)
calculate_ubl(e5ds,ID=ID,days=1)
calculate_ubl(e5ds,ID=ID,days=3)
calculate_ubl(e5ds,ID=ID,days=7)
calculate_ubl(e5ds,ID=ID,days=30)

calculate_vtbl(e5ds,ID=ID)
calculate_vtbl(e5ds,ID=ID,days=1)
calculate_vtbl(e5ds,ID=ID,days=3)
calculate_vtbl(e5ds,ID=ID,days=7)
calculate_vtbl(e5ds,ID=ID,days=30)

calculate_vbl(e5ds,ID=ID)
calculate_vbl(e5ds,ID=ID,days=1)
calculate_vbl(e5ds,ID=ID,days=3)
calculate_vbl(e5ds,ID=ID,days=7)
calculate_vbl(e5ds,ID=ID,days=30)