using DrWatson
@quickactivate "S2DExploration"
using ERA5Reanalysis
using Statistics

function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

function smoothing!(
    data :: Vector{FT};
    days :: Int
) where FT <: Real

    ndt = length(data)
    buffer,weights = calculatebufferweights(days*24)

    for idt = 1 : (ndt-2*buffer)
        iidata = @views data[idt .+ (0:(2*buffer))]
        data[idt] = sum(iidata.*weights)
    end

    for idt = 1 : (ndt-2*buffer)
        data[idt+buffer] = data[idt,ilvl]
    end

    data[1:buffer] .= NaN
    data[(ndt-buffer+1):end] .= NaN

    return

end

function smoothing!(
    data :: Array{FT,2};
    days :: Int
) where FT <: Real

    ndt,nlvl = size(data)
    buffer,weights = calculatebufferweights(days*24)

    for ilvl = 1 : nlvl, idt = 1 : (ndt-2*buffer)
        iidata = @views data[idt .+ (0:(2*buffer)),ilvl]
        data[idt,ilvl] = sum(iidata.*weights)
    end

    for ilvl = 1 : nlvl, idt = 1 : (ndt-2*buffer)
        data[idt+buffer,ilvl] = data[idt,ilvl]
    end

    data[1:buffer,:] .= NaN
    data[(ndt-buffer+1):end,:] .= NaN

    return

end