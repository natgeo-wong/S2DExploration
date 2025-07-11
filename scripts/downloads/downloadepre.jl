using DrWatson
@quickactivate "S2DExploration"
using ERA5Reanalysis

e5ds = ERA5Hourly(start=Date(1991,1,1),stop=Date(2024,12,31),path=datadir())
ereg = ERA5Region("KS_OK",path=srcdir())
ereg = ERA5Region("TN_MS_AL",path=srcdir())

epre = PressureVariable("u",hPa=0); download(e5ds,epre,ereg,pall=true)
epre = PressureVariable("v",hPa=0); download(e5ds,epre,ereg,pall=true)
epre = PressureVariable("w",hPa=0); download(e5ds,epre,ereg,pall=true)
epre = PressureVariable("t",hPa=0); download(e5ds,epre,ereg,pall=true)
epre = PressureVariable("q",hPa=0); download(e5ds,epre,ereg,pall=true)
epre = PressureVariable("cc",hPa=0); download(e5ds,epre,ereg,pall=true)
