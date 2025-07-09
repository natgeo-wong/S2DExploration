using DrWatson
@quickactivate "S2DExploration"
using ERA5Reanalysis

e5ds = ERA5Hourly(start=Date(1991,1,1),stop=Date(2024,12,31),path=datadir())
ereg = ERA5Region("SGP",path=srcdir())
ereg = ERA5Region("BNF",path=srcdir())

esgl = [
    SingleVariable("tsr"),SingleVariable("ttr"),SingleVariable("ssr"),
    SingleVariable("slhf"),SingleVariable("sshf"),SingleVariable("str"),
    SingleVariable("sst"),SingleVariable("t2m"),SingleVariable("skt"),
    SingleVariable("tcw"),SingleVariable("tcwv"),SingleVariable("sp"),
    SingleVariable("hcc"),SingleVariable("mcc"),SingleVariable("lcc"),
    SingleVariable("tcc"),
]
download(e5ds,esgl,ereg)