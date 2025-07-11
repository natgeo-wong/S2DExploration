using DrWatson
@quickactivate "S2DExploration"

using GeoRegions
using GeoInterface
using JSON3

json_string = read(datadir("us-state-boundaries.json"),String);
data = JSON3.read(json_string);

## Retrieve Oklahoma State Boundary Coordinates
coords = data[11].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end

geo1 = GeoRegion(lon,lat,ID = "OK", pID = "GLB", name = "Oklahoma")

## Retrieve Kansas State Boundary Coordinates
coords = data[20].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end
geo2 = GeoRegion(lon,lat,ID = "KS", pID = "GLB", name = "Kansas")

## Retrieve Alabama State Boundary Coordinates
coords = data[25].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end
geo3 = GeoRegion(lon,lat,ID = "AL", pID = "GLB", name = "Alabama")

## Retrieve Georgia State Boundary Coordinates
coords = data[5].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end
geo4 = GeoRegion(lon,lat,ID = "GA", pID = "GLB", name = "Georgia")

## Retrieve Tennessee State Boundary Coordinates
coords = data[47].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end
geo5 = GeoRegion(lon,lat,ID = "TN", pID = "GLB", name = "Tennessee")

## Retrieve Mississippi State Boundary Coordinates
coords = data[52].st_asgeojson.geometry.coordinates[1]; ncoord = length(coords)
lon = zeros(ncoord)
lat = zeros(ncoord)

for ii = 1 : ncoord
    lon[ii] = coords[ii][1]
    lat[ii] = coords[ii][2]
end
geo6 = GeoRegion(lon,lat,ID = "MS", pID = "GLB", name = "Mississippi")

a = union(geo1.geometry.polygon,geo2.geometry.polygon)
# b = union(geo3.geometry.polygon,geo4.geometry.polygon,geo5.geometry.polygon,geo6.geometry.polygon)

lg1 = GeoInterface.convert(LibGEOS, geo1.geometry.polygon)
lg2 = GeoInterface.convert(LibGEOS, geo2.geometry.polygon)
a = LibGEOS.union(lg1, lg2)
a = GeoInterface.convert(GeometryBasics, a)
apnts = GeoInterface.coordinates(a)[1]; npnts = length(apnts)
lon = zeros(npnts)
lat = zeros(npnts)

for ipnt = 1 : npnts
    lon[ipnt] = apnts[ipnt][1]
    lat[ipnt] = apnts[ipnt][2]
end

rmID("KS_OK",path=srcdir())
geoa = GeoRegion(lon,lat,ID = "KS_OK",pID="GLB",name="Test",save = true,path=srcdir())

lg3 = GeoInterface.convert(LibGEOS, geo3.geometry.polygon)
lg4 = GeoInterface.convert(LibGEOS, geo4.geometry.polygon)
lg5 = GeoInterface.convert(LibGEOS, geo5.geometry.polygon)
lg6 = GeoInterface.convert(LibGEOS, geo6.geometry.polygon)
b = LibGEOS.union(lg3,lg5)
b = LibGEOS.union(b,lg6)
b = GeoInterface.convert(GeometryBasics, b)
apnts = GeoInterface.coordinates(b)[1]; npnts = length(apnts)
lon = zeros(npnts)
lat = zeros(npnts)

for ipnt = 1 : npnts
    lon[ipnt] = apnts[ipnt][1]
    lat[ipnt] = apnts[ipnt][2]
end

rmID("TN_MS_AL",path=srcdir())
geoa = GeoRegion(lon,lat,ID = "TN_MS_AL",pID="GLB",name="Test",save = true,path=srcdir())