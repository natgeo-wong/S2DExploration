using DrWatson
@quickactivate "S2DExploration"
using NASAPrecipitation

npd = IMERGFinalHH(start=Date(2001),stop=Date(2001),path=datadir())
geo = GeoRegion("KS_OK",path=srcdir())

download(npd,geo)