function stats_df(df, ts, f, dt)
    nam = names(df)
    mapfoldl(vcat, ts) do tts
        @chain begin
            filter(:t => t -> (t .>= tts) .& (t .< tts + Minute(dt)), df)
            map(f, eachcol(_[!, 2:end]))
            hcat(DataFrame(; t = tts), DataFrame(_', nam[2:end]))
        end
    end
end

function val2zeroD(x::T) where {T}
    y = Array{T}(undef)
    return fill!(y, x)
end

function nan2zero(x)
    y = deepcopy(x)
    x[isnan.(x)] .= 0.0
    return x
end

function get_DMA_config(Qsh::Float64, Qsa::Float64, T::Float64, p::Float64, column::Symbol)
    lpm = 1.666666e-5
    (column == :TSILONG) && ((r₁, r₂, l) = (9.37e-3, 1.961e-2, 0.44369))
    (column == :HFDMA) && ((r₁, r₂, l) = (0.05, 0.058, 0.6))
    (column == :RDMA) && ((r₁, r₂, l) = (2.4e-3, 50.4e-3, 10e-3))
    (column == :HELSINKI) && ((r₁, r₂, l) = (2.65e-2, 3.3e-2, 10.9e-2))
    (column == :VIENNASHORT) && ((r₁, r₂, l) = (25e-3, 33.5e-2, 0.11))
    (column == :VIENNAMEDIUM) && ((r₁, r₂, l) = (25e-3, 33.5e-2, 0.28))
    (column == :VIENNALONG) && ((r₁, r₂, l) = (25e-3, 33.5e-2, 0.50))

    form = (column == :RDMA) ? :radial : :cylindrical

    qsh = Qsh * lpm
    qsa = Qsa * lpm
    t = T + 273.15
    leff = 13.0
    polarity = :+

    Λ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, polarity, 6, form)

    return Λ
end

vtod(Λ, v) = @chain vtoz(Λ, v) ztod(Λ, 1, _)

function regrid(N, D, δ)
    return map(1:length(δ.Dp)) do i
        ii = (D .< δ.De[i]) .& (D .> δ.De[i+1])
        return mean(N[ii])
    end
end

getf64(x, i) = @match x begin
    "missing" => missing
    _ => @chain split(x, ",") _[i] parse.(Float64, _)
end

getRIE(x, i) = @match x begin
    "missing" => missing
    _ => @chain split(x, ",") _[4] parse(Int16, _, base = 16)
end

function readraw(ss)
    bp = "/home/mdpetter/epcape/mtsncsudenudedccn/"

    file = @chain begin
        files = readdir(bp)
        filter(f -> split(f, ".")[end] == "txt", _)
        filter(_) do f
            a = split.(f, "_")[2]
            d = Date(split(a, ".")[1], dateformat"yyyymmdd")
            return (d >= ss) & (d <= ss)
        end
    end

    s = open(bp * file[1]) do file
        return read(file, String)
    end

    a = split(s, "\n")
    head = map(s -> split(s, ";")[1], a[1:end-1])
    tcpc = map(s -> split(s, ";")[2], a[1:end-1])
    tpops = map(s -> split(s, ";")[3], a[1:end-1])
    tccn = map(s -> split(s, ";")[4], a[1:end-1])
    tlj = map(s -> split(s, ";")[5], a[1:end-1])

    # General
    t = map(
        x ->
            (@chain split(x, ";") getindex(_, 1) split(_, ",") getindex(_, 1) DateTime),
        head,
    )
    denuded = map(
        x -> (@chain split(x, ";") getindex(_, 1) split(_, ",") getindex(_, 2) parse(
            Bool,
            _,
        )),
        head,
    )

    return t, denuded, tcpc, tpops, tccn, tlj
end

function time2cdf(cdfname, t, ss)
    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))
    base_time = Array{Int}(undef)
    fill!(base_time, datetime2unix(DateTime(ss)))
    ncwrite(base_time, cdfname, "base_time")
    offset = @chain collect(ts) datetime2unix.(_) Int32.(_) .-(_, base_time)
    ncwrite(offset, cdfname, "time_offset")

    return nothing
end

function pops2cdf(cdfname, t, tpops, ss)
    DminPOPS = [
        158,
        175,
        194,
        215,
        238,
        292,
        381,
        500,
        592,
        820,
       1290,
       1674,
       2090,
       2917,
       3795,
       4923,
    ]
    DmaxPOPS = [
        175,
        194,
        215,
        238,
        292,
        381,
        500,
        592,
        820,
       1290,
       1674,
       2090,
       2917,
       3795,
       4923,
       6393,
    ]
    DboundPOPS = hcat(DminPOPS, DmaxPOPS)'[:, :]
    DmidPOPS = DminPOPS .+ (DmaxPOPS .- DminPOPS) ./ 2.0
    dlogDpPOPS = log.(DmaxPOPS ./ DminPOPS)

    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))

    qpops = map(x -> getf64(x, 8), tpops)
    popscounts = map(x -> getf64(x, 12:27) ./ dlogDpPOPS, tpops)
    ii = .~ismissing.(qpops)
    popspsd = DataFrame(hcat(popscounts[ii]...)'[:, :] ./ qpops[ii], :auto)
    popsn = map(x -> sum(Vector(x) .* dlogDpPOPS), eachrow(popspsd))
    popsdf = @chain begin
        hcat(DataFrame(; t = t[ii], Nt = popsn, q = qpops[ii]), popspsd)
        stats_df(_, ts, mean, dtav)
    end

    ncwrite(popsdf[!, :Nt], cdfname, "total_N_conc_pops")
    ncwrite(popsdf[!, :q], cdfname, "sample_flow_pops")
    ncwrite(Matrix(popsdf[!, 4:end])'[:, :], cdfname, "dN_dlog_Dp_POPS")
    ncwrite(DmidPOPS, cdfname, "diameter_pops_midpoints")
    ncwrite(DboundPOPS, cdfname, "diameter_pops_bounds")

    return nothing
end

function state2cdf(cdfname, t, tlj, ss)
    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))
    RHdma = map(x -> getf64(x, 7), tlj)
    Tdma = map(x -> getf64(x, 8), tlj)
    statedf = @chain begin
        DataFrame(; t = t, T = Tdma, RH = RHdma)
        stats_df(_, ts, mean, dtav)
    end

    ncwrite(statedf[!, :T], cdfname, "sample_temperature")
    ncwrite(statedf[!, :RH], cdfname, "sample_relative_humidity")

    return nothing
end

function dmasetup2cdf(cdfname, Λ, δ)
    Dbound = hcat(δ.De[2:end], δ.De[1:end-1])
    x = (Λ.polarity == :+) ? Int32(1) : Int32(0)
    ncwrite(δ.𝐀[:, :], cdfname, "inversion_matrix")
    ncwrite(δ.Dp, cdfname, "diameter_mobility_midpoints")
    ncwrite(Dbound'[:, :], cdfname, "diameter_mobility_bounds")
    ncwrite(val2zeroD(Λ.r1), cdfname, "DMA_inner_radius")
    ncwrite(val2zeroD(Λ.r2), cdfname, "DMA_outer_radius")
    ncwrite(val2zeroD(Λ.l), cdfname, "DMA_length")
    ncwrite(val2zeroD(Λ.leff), cdfname, "DMA_effective_length")
    ncwrite(val2zeroD(x), cdfname, "hv_polarity")
    ncwrite(val2zeroD(Λ.qsh), cdfname, "sheath_flow_dma")
    ncwrite(val2zeroD(Λ.qsa), cdfname, "sample_flow_dma")
    ncwrite(val2zeroD(Λ.t), cdfname, "reference_gas_temperature")
    ncwrite(val2zeroD(Λ.p), cdfname, "reference_gas_pressure")

    return nothing
end

function ccnflow2cdf(cdfname, t, tccn, ss)
    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))
    qsa = map(x -> getf64(x, 12), tccn)
    qsh = map(x -> getf64(x, 13), tccn)

    ccndf = @chain begin
        DataFrame(; t = t, qsh = qsh, qsa = qsa)
        stats_df(_, ts, mean, dtav)
        coalesce.(_, NaN)
    end

    ncwrite(ccndf[!, :qsa], cdfname, "sample_flow_ccn")
    ncwrite(ccndf[!, :qsh], cdfname, "sheath_flow_ccn")

    return nothing
end

function latlonalt2cdf(cdfname)
    ncwrite(val2zeroD(32.8405), cdfname, "lat")
    ncwrite(val2zeroD(-117.2496), cdfname, "lon")
    return ncwrite(val2zeroD(220.0), cdfname, "alt")
end

function distribution2cdf(cdfname, t, tcpc, tccn, tlj, denuded, Λ, δ, ss)
    cpcflag = map(x -> getRIE(x, 4), tcpc)
    ccn = map(x -> getf64(x, 21), tccn)
    dsd = map(x -> getf64(x, 23:42), tccn)
    T1 = map(x -> getf64(x, 5), tccn)
    T2 = map(x -> getf64(x, 6), tccn)
    T3 = map(x -> getf64(x, 7), tccn)
    dT = map(x -> getf64(x, 17), tccn)

    # LJ decode    
    state = map(x -> (@chain split(x, ",") getindex(_, 4)), tlj)
    Vr = map(x -> getf64(x, 5), tlj)
    cvi_raw = map(x -> getf64(x, 6), tlj)
    cvi = map(x -> (x > -20.0) ? 0.0 : 1.0, cvi_raw)
    if ss < Date(2023,3,10)
        cvi[:] .= 0.0
    end
    ccount = map(x -> getf64(x, 9), tlj)
    cserial = map(x -> getf64(x, 3), tcpc)

    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))

    rcn_serial = NaN .* ones(288, 100)
    rcn_count = NaN .* ones(288, 100)
    rccn = NaN .* ones(288, 100)
    invcnserial_0 = NaN .* ones(288, 100)
    invcncount_0 = NaN .* ones(288, 100)
    invccn_0 = NaN .* ones(288, 100)
    invcnserial_2 = NaN .* ones(288, 100)
    invcncount_2 = NaN .* ones(288, 100)
    invccn_2 = NaN .* ones(288, 100)
    meandsd = NaN .* ones(288, 20)
    ccn_dT = NaN .* ones(288)
    isdenuded = NaN .* falses(288)
    iscvi = NaN .* ones(288)
    cpc_status = NaN .* ones(Int, 288)
    scans = [1; findall(denuded[1:end-1] .⊻ denuded[2:end])]

    i = 1
    map(1:length(scans)-1) do i
        ii = scans[i]+4:scans[i+1]
        n = length(ii)
        if n < 590
            println("Scan Excluded")
            return nothing
        end

        jj = state[ii] .== "UPSCAN"
        kk = state[ii] .== "DOWNSCAN"
        tx = round((t[ii][1]), Dates.Minute)
        ti = findall(tx .== ts)[1]

        Djj = map(v -> vtod(Λ, v), Vr[ii][jj])
        Dkk = map(v -> vtod(Λ, v), Vr[ii][kk])

        isdenuded[ti] = denuded[ii][10]
        isdenuded[ti+1] = denuded[ii][10]
        iscvi[ti] = mean(cvi[ii][jj])
        iscvi[ti+1] = mean(cvi[ii][kk])
        ccn_dT[ti] = dT[ii][10]
        ccn_dT[ti+1] = dT[ii][10]
        cpc_stats = cpcflag[ii][10]

        delay_time_cpc_serial = 5
        delay_time_cpc_count = 4
        delay_time_ccn = (denuded[ii][10] == true) ? 41 : 21

        cserial_shift = circshift(cserial, -delay_time_cpc_serial)
        ccount_shift = circshift(ccount, -delay_time_cpc_count)
        ccn_shift = circshift(ccn, -delay_time_ccn)
        dsd_shift = circshift(ccn, -delay_time_ccn)

        cnserialjj = @chain cserial_shift[ii][jj] regrid(_, Djj, δ)
        cncountjj = @chain ccount_shift[ii][jj] regrid(_, Djj, δ)
        ccnjj = @chain ccn_shift[ii][jj] regrid(_, Djj, δ)
        cnserialkk = @chain cserial_shift[ii][kk] regrid(_, Dkk, δ)
        cncountkk = @chain ccount_shift[ii][kk] regrid(_, Dkk, δ)
        ccnkk = @chain ccn_shift[ii][kk] regrid(_, Dkk, δ)
        meandsdjj = mean(hcat(dsd[ii][jj]...); dims = 2)[:]
        meandsdkk = mean(hcat(dsd[ii][kk]...); dims = 2)[:]

        rcn_serial[ti, :] = cnserialjj
        rcn_serial[ti+1, :] = cnserialkk
        rcn_count[ti, :] = cncountjj
        rcn_count[ti+1, :] = cncountkk
        rccn[ti, :] = ccnjj
        rccn[ti+1, :] = ccnkk
        meandsd[ti, :] = meandsdjj
        meandsd[ti+1, :] = meandsdkk

        invcnserialjj_2 =
            rinv2(nan2zero(cnserialjj), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)
        invcnserialjj_0 =
            rinv2(nan2zero(cnserialjj), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invcncountjj_2 =
            rinv2(nan2zero(cncountjj), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)
        invcncountjj_0 =
            rinv2(nan2zero(cncountjj), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invccnjj_0 =
            rinv2(nan2zero(ccnjj), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invccnjj_2 =
            rinv2(nan2zero(ccnjj), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)

        invcnserialkk_2 =
            rinv2(nan2zero(cnserialkk), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)
        invcnserialkk_0 =
            rinv2(nan2zero(cnserialkk), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invcncountkk_2 =
            rinv2(nan2zero(cncountkk), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)
        invcncountkk_0 =
            rinv2(nan2zero(cncountkk), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invccnkk_0 =
            rinv2(nan2zero(ccnkk), δ; λ₁ = 0.1, λ₂ = 20.0, order = 0, initial = true)
        invccnkk_2 =
            rinv2(nan2zero(ccnkk), δ; λ₁ = 5.0, λ₂ = 20.0, order = 2, initial = false)

        invcnserial_0[ti, :] = invcnserialjj_0.S
        invcnserial_0[ti+1, :] = invcnserialkk_0.S
        invcncount_0[ti, :] = invcncountjj_0.S
        invcncount_0[ti+1, :] = invcncountkk_0.S
        invccn_0[ti, :] = invccnjj_0.S
        invccn_0[ti+1, :] = invccnkk_0.S

        invcnserial_2[ti, :] = invcnserialjj_2.S
        invcnserial_2[ti+1, :] = invcnserialkk_2.S
        invcncount_2[ti, :] = invcncountjj_2.S
        invcncount_2[ti+1, :] = invcncountkk_2.S
        invccn_2[ti, :] = invccnjj_2.S
        invccn_2[ti+1, :] = invccnkk_2.S
        return println(tx)
    end

    ccn_ss = round.(0.0777 .* ccn_dT .- 0.14895; digits = 1)

    ncwrite(Int32.(isdenuded), cdfname, "denuder_flag")
    ncwrite(iscvi, cdfname, "cvi_flag")
    ncwrite(ccn_dT, cdfname, "ccn_gradient")
    ncwrite(ccn_ss, cdfname, "ccn_supersaturation")
    ncwrite(nan2zero(cpc_status), cdfname, "qc_cpc_flag")
    ncwrite(rcn_serial'[:, :], cdfname, "response_function_cpc_serial")
    ncwrite(rcn_count'[:, :], cdfname, "response_function_cpc_count")
    ncwrite(rccn'[:, :], cdfname, "response_function_ccn")
    ncwrite(invcnserial_0'[:, :], cdfname, "dN_dlog_Dp_DMA_serial_0")
    ncwrite(invcnserial_2'[:, :], cdfname, "dN_dlog_Dp_DMA_serial_2")
    ncwrite(invcncount_0'[:, :], cdfname, "dN_dlog_Dp_DMA_count_0")
    ncwrite(invcncount_2'[:, :], cdfname, "dN_dlog_Dp_DMA_count_2")
    ncwrite(invccn_0'[:, :], cdfname, "dN_dlog_Dp_CCN_0")
    ncwrite(invccn_2'[:, :], cdfname, "dN_dlog_Dp_CCN_2")
    ncwrite(meandsd'[:, :], cdfname, "droplet_size_distribution_ccn")

    return nothing
end

function fitme(
    ccn_supersaturation,
    diameter_mobility_midpoints,
    response_function_ccn,
    response_function_cpc_count,
    i,
)
    ig(x) = @match x begin
        0.2 => 100.0
        0.4 => 58.0
        0.6 => 47.0
        0.8 => 40.0
        1.0 => 35.0
        _ => 50.0
    end

    Taf(Dp, μ, σ) = @. 0.5 * (1.0 + erf((Dp - μ) ./ (sqrt(2.0 .* σ))))
    model(x, p) = p[3] * Taf(x, p[1], p[2])

    function threshold(N, Dp, c::Float64, n1::Float64, n2::Float64)
        M = deepcopy(N)
        M[(N.<=c).&(Dp.>200)] .= n2
        M[(N.<=c).&(Dp.<200)] .= n1
        return M
    end

    mig = ig.(ccn_supersaturation)

    initial = ([mig[i], mig[i] ./ 10.0, 0.9])

    tCCN =
        threshold(response_function_ccn[:, i], diameter_mobility_midpoints, 0.1, 0.1, 0.1)
    tCN = threshold(
        response_function_cpc_count[:, i],
        diameter_mobility_midpoints,
        0.1,
        0.0,
        0.1,
    )
    af = tCCN ./ tCN
    jj = (diameter_mobility_midpoints .< 180.0) .& isfinite.(af)

    Ax = try
        fit = curve_fit(model, diameter_mobility_midpoints[jj], af[jj], initial)
        if fit.converged == true
            fit.param
        else
            [1000.0, 100.0, 1.0]
        end
    catch
         [1000.0, 100.0, 1.0]
    end

    p = scatter(
        diameter_mobility_midpoints[jj],
        af[jj];
        ylim = (0, 1.4),
        xscale = :log10,
        xlim = (10, 400),
        minorgrid = true,
        label = "rCCN/rCN $(ccn_supersaturation[i])",
        legend = :topleft,
    )

    κ = map(ccn_supersaturation[i], Ax[1] * 1e-9) do sc, Dd
        k = try
            get_κ(sc, Dd)
        catch
            NaN
        end
        o = (Dd == 1e-6) ? NaN : k
        o
    end

    return p = plot!(
        diameter_mobility_midpoints,
        model(diameter_mobility_midpoints, Ax);
        label = "$(round(Ax[1]))",
    ),
    κ
end
