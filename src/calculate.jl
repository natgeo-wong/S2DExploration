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

Tb2OLR(Tb::Real) = 5.67e-8 * (Tb * (1.228 - 1.106e-3*Tb))^4

function calculate_ptbl(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar = PressureVariable("z")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    z = ds[evar.ncID][:,:] ./ 9.80665; ndt = size(z,2)
    close(ds)

    evar = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("blh")
    ds = read_climatology(ID,e5ds,evar,days=days)
    blh = ds[evar.ncID][:]
    close(ds)
    
    ptbl = zeros(ndt)

    for idt = 1 : ndt

        isp = sp[idt]
        ip = vcat(p[(p.<isp)],isp)
        iz = vcat(z[(p.<isp),idt],0)

        spl = Spline1D(reverse(iz),reverse(ip),k=1)
        iptbl = spl(blh[idt])
        ptbl[idt] = iptbl > isp ? isp : iptbl

    end

    save_climatology(ID,e5ds,SingleVariable("ptbl",path=srcdir()),ptbl,days=days)

end

function calculate_utbl(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar = PressureVariable("u")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    u = ds[evar.ncID][:,:]; ndt = size(u,2)
    close(ds)

    evar = SingleVariable("u10")
    ds = read_climatology(ID,e5ds,evar,days=days)
    u10 = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("ptbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:]
    close(ds)
    
    utbl = zeros(ndt)

    for idt = 1 : ndt

        isp = sp[idt]
        ip = vcat(p[(p.<isp)],    isp)
        iu = vcat(u[(p.<isp),idt],u10[idt])

        spl = Spline1D(ip,iu,k=1)
        utbl[idt] = spl(ptbl[idt])

    end

    save_climatology(ID,e5ds,SingleVariable("utbl",path=srcdir()),utbl,days=days)


end

function calculate_vtbl(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar = PressureVariable("v")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    v = ds[evar.ncID][:,:]; ndt = size(v,2)
    close(ds)

    evar = SingleVariable("v10")
    ds = read_climatology(ID,e5ds,evar,days=days)
    v10 = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("ptbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:]
    close(ds)
    
    vtbl = zeros(ndt)

    for idt = 1 : ndt

        isp = sp[idt]
        ip = vcat(p[(p.<isp)],    isp)
        iv = vcat(v[(p.<isp),idt],v10[idt])

        spl = Spline1D(ip,iv,k=1)
        vtbl[idt] = spl(ptbl[idt])

    end

    save_climatology(ID,e5ds,SingleVariable("vtbl",path=srcdir()),vtbl,days=days)


end

function calculate_ubl(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar = PressureVariable("u")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    u = ds[evar.ncID][:,:]; ndt = size(u,2)
    close(ds)

    evar = SingleVariable("u10")
    ds = read_climatology(ID,e5ds,evar,days=days)
    u10 = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("utbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    utbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("ptbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:]
    close(ds)
    
    ubl = zeros(ndt)

    for idt = 1 : ndt

        isp = sp[idt]
        iptbl = ptbl[idt]
        ip = vcat(ptbl[idt],p[(p.<isp).&(p.>iptbl)],    isp); np = length(ip)
        iu = vcat(utbl[idt],u[(p.<isp).&(p.>iptbl),idt],u10[idt])

        ubl[idt] = trapz(ip,iu) / trapz(ip,ones(np))

    end

    save_climatology(ID,e5ds,SingleVariable("ubl",path=srcdir()),ubl,days=days)


end

function calculate_vbl(
    e5ds  :: ERA5Hourly;
    ID    :: String,
    days  :: Int = 0
)

    evar = PressureVariable("v")
    ds = read_climatology(ID,e5ds,evar,days=days)
    p = ds["pressures"][:] * 100
    v = ds[evar.ncID][:,:]; ndt = size(v,2)
    close(ds)

    evar = SingleVariable("v10")
    ds = read_climatology(ID,e5ds,evar,days=days)
    v10 = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("vtbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    vtbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("ptbl",path=srcdir())
    ds = read_climatology(ID,e5ds,evar,days=days)
    ptbl = ds[evar.ncID][:]
    close(ds)

    evar = SingleVariable("sp")
    ds = read_climatology(ID,e5ds,evar,days=days)
    sp = ds[evar.ncID][:]
    close(ds)
    
    vbl = zeros(ndt)

    for idt = 1 : ndt

        isp = sp[idt]
        iptbl = ptbl[idt]
        ip = vcat(ptbl[idt],p[(p.<isp).&(p.>iptbl)],    isp); np = length(ip)
        iv = vcat(vtbl[idt],v[(p.<isp).&(p.>iptbl),idt],v10[idt])

        vbl[idt] = trapz(ip,iv) / trapz(ip,ones(np))

    end

    save_climatology(ID,e5ds,SingleVariable("vbl",path=srcdir()),vbl,days=days)


end