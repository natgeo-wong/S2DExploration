using DrWatson
@quickactivate "S2DExploration"
using NASAMergedTb

tbd = TbDataset(start=Date(2001),stop=Date(2001),path=datadir())
geo = GeoRegion("KS_OK",path=srcdir())

download(tbd,geo,overwrite=true)