using DrWatson
using Dierckx
using ERA5Reanalysis
using Statistics
using Trapz

include(srcdir("common.jl"))

function calculate_ptrop(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar  = PressureVariable("t")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:]; iip = (p.>=20) .& (p .<= 200)
    p = p[iip] * 100
    t = ds[evar.ncID][iip,:]; ndt = size(t,2)
    close(ds)

    pp = 2000 : 100 : 20000; np = length(pp)
    dtdp = zeros(np); ptrop = zeros(ndt)

    for idt = 1 : ndt

        iit = @views t[:,idt]
        spl = Spline1D(p,iit)

        for ip in 1 : np
            iipp = pp[ip]
            dtdp[ip] = derivative(spl,iipp) * 9.81 * iipp / 287 / evaluate(spl,iipp)
        end

        ip = np; idtdp = dtdp[ip]
        while idtdp > 0.002
            ip -= 1
            idtdp = dtdp[ip]
        end
        ptrop[idt] = pp[ip]

    end

    evar = SingleVariable("ptrop",path=srcdir())
    save_climatology(ID,e5ds,evar,ptrop,days=days)
    
end

function calculate_ωc(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar  = PressureVariable("w")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:]
    w = ds[evar.ncID][:,:]; ndt = size(t,2)
    close(ds)

    evar  = PressureVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:] / 100
    close(ds)

    evar  = PressureVariable("ptrop")
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptrop = ds[evar.ncID][:] / 100
    close(ds)

    nωc = 2; ωc = zeros(nωc+1,idt)


    for idt = 1 : ndt

        ip = @views p[(p.>ptrop[idt]).&(p.<sp[idt])]
        iw = @views w[(p.>ptrop[idt]).&(p.<sp[idt]),idt]
        ix = vcat(ip,sp)
        iy = vcat(iw,0)
        spl = Spline1D(ix,iy)
        wspl = zeros(101)
        pspl = zeros(101)
        pmin = minimum(ip)

        for ii = 1 : 100

            iip = pmin + (sp[idt]-pmin) * (ii-1) / 99
            pspl[ii+1] = iip
            wspl[ii+1] = evaluate(spl,iip)

        end
        pspl[1] = ptrop[idt]

        for ic = 1 : nωc
            ωc[ic,idt] = trapz(
                pspl .- pspl[1],
                wspl .* sin.(ic*pi*(pspl .- pspl[1]) / (pspl[end] - pspl[1]))
            ) * 2 / (pspl[end] - pspl[1])
        end

    end

    evar = SingleVariable("ω1",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[1,:])
    evar = SingleVariable("ω2",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[2,:])
    # evar = SingleVariable("ωs",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[3,:])
    
end