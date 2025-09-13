using DrWatson
using Dierckx
using ERA5Reanalysis
using Logging
using Statistics
using Trapz

include(srcdir("common.jl"))

function calculate_dtdp(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar  = PressureVariable("t")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100; np = length(p)
    t = ds[evar.ncID][:,:]; ndt = size(t,2)
    close(ds)

    dtdp = zeros(np,ndt)

    for idt = 1 : ndt

        iit = @views t[:,idt]
        spl = Spline1D(p,iit)

        for ip in 1 : np
            iipp = p[ip]
            dtdp[ip,idt] = derivative(spl,iipp) * 9.81 * iipp / 287 / evaluate(spl,iipp)
        end

    end

    evar = PressureVariable("dtdp",path=srcdir())
    save_climatology(ID,e5ds,evar,dtdp,Int.(p./100),days=days)
    
end

function calculate_dtdp_hres(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar  = PressureVariable("t")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    t = ds[evar.ncID][:,:]; ndt = size(t,2)
    close(ds)

    pp = 1000 : 1000 : 100000; np = length(pp)

    dtdp = zeros(np,ndt)

    for idt = 1 : ndt

        iit = @views t[:,idt]
        spl = Spline1D(p,iit)

        for ip in 1 : np
            iipp = pp[ip]
            dtdp[ip,idt] = derivative(spl,iipp) * 9.81 * iipp / 287 / evaluate(spl,iipp)
        end

    end

    evar = PressureVariable("dtdp_hres",path=srcdir())
    save_climatology(ID,e5ds,evar,dtdp,Int.(pp./100),days=days)
    
end

function calculate_ptrop(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar  = PressureVariable("t")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    t = ds[evar.ncID][:,:]; ndt = size(t,2)
    close(ds)

    pp = 2000 : 1000 : 50000; np = length(pp)
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
    w = ds[evar.ncID][:,:]; ndt = size(w,2)
    close(ds)

    evar  = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:] / 100
    close(ds)

    evar  = SingleVariable("ptrop",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptrop = ds[evar.ncID][:] / 100
    close(ds)

    nωc = 2; ωc = zeros(nωc+1,ndt)

    for idt = 1 : ndt

        if !isnan(sp[idt])

            ip = @views p[(p.>ptrop[idt]).&(p.<sp[idt])]
            iw = @views w[(p.>ptrop[idt]).&(p.<sp[idt]),idt]
            ix = vcat(ip,sp[idt])
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

    end

    evar = SingleVariable("ω1",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[1,:],days=days)
    evar = SingleVariable("ω2",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[2,:],days=days)
    # evar = SingleVariable("ωs",path=srcdir()); save_climatology(ID,e5ds,evar,ωc[3,:],days=days)
    
end