using Plots
using Plots.PlotMeasures
using LsqFit
using MLStyle
using SpecialFunctions
using Chain

include("reader.jl")
include("rh_functions.jl")
include("processingroutines.jl")

function quicklook(week_to_analyze)
    tf,
    total_N_conc_pops,
    sample_flow_pops,
    dNdlog_Dp_POPS,
    diameter_pops_midpoints,
    diameter_pops_bounds,
    sample_temperature,
    sample_relative_humidity,
    inversion_matrix,
    diameter_mobility_midpoints,
    diameter_mobility_bounds,
    qc_cpc_flag,
    ccn_supersaturation,
    denuder_flag,
    response_function_cpc_serial,
    response_function_cpc_count,
    response_function_ccn,
    dN_dlog_Dp_DMA_serial_0,
    dN_dlog_Dp_DMA_serial_2,
    dN_dlog_Dp_DMA_count_0,
    dN_dlog_Dp_DMA_count_2,
    dN_dlog_Dp_CCN_0,
    dN_dlog_Dp_CCN_2,
    droplet_size_distribution_ccn,
    sample_flow_ccn,
    sheath_flow_ccn,
    cvi_flag =
        read_netCDF("/home/mdpetter/epcape/epc10denudedccnS2.b1/", week_to_analyze)

    ss = tf[1]
    ee = tf[1] + Day(7)
    tm = ss:Hour(24):ee
    ticks = Dates.format.(tm, "mm/dd")

    p1 = plot(
        tf,
        sample_temperature;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylabel = "T (°C)",
    )

    p2 = plot(
        tf,
        sample_relative_humidity;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylabel = "RH (%)",
    )

    p3 = plot(
        tf,
        sample_flow_pops;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylabel = "Qpops (ccm)",
    )

    p4 = plot(
        tf,
        sheath_flow_ccn .+ sample_flow_ccn;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylabel = "Q CCN (ml min-1)",
    )

    p5 = plot(
        tf,
        sheath_flow_ccn ./ sample_flow_ccn;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylabel = "Q sh-to-sa (-)",
    )

    p6 = plot(
        tf,
        cvi_flag;
        xticks = (tm, ticks),
        ylim = (0,1.2),
        minorgrid = true,
        color = :black,
        label = :none,
        ylabel = "CVI Flag (-)",
    )

    phk =
        plot(p1, p2, p3, p4, p5, p6; layout = grid(6, 1), xlim = (ss, ee), size = (900, 800), left_margin = 60px, right_marking = 40px)

    dstr = Dates.format(week_to_analyze, "yyyymmdd")
    savefig(phk, "src/assets/hk/" * dstr * "_hk.png")

    dlnDp = log.(diameter_mobility_bounds[2, :] ./ diameter_mobility_bounds[1, :])

    N_count0 = map(i -> sum(dN_dlog_Dp_DMA_count_0[:, i] .* dlnDp), 1:length(tf))
    N_count2 = map(i -> sum(dN_dlog_Dp_DMA_count_2[:, i] .* dlnDp), 1:length(tf))
    N_serial0 = map(i -> sum(dN_dlog_Dp_DMA_serial_0[:, i] .* dlnDp), 1:length(tf))
    N_serial2 = map(i -> sum(dN_dlog_Dp_DMA_serial_2[:, i] .* dlnDp), 1:length(tf))

    N_ccn0 = map(i -> sum(dN_dlog_Dp_CCN_0[:, i] .* dlnDp), 1:length(tf))
    N_ccn2 = map(i -> sum(dN_dlog_Dp_CCN_2[:, i] .* dlnDp), 1:length(tf))

    kk = diameter_mobility_midpoints .> 130.0

    p1 = plot(
        tf,
        total_N_conc_pops;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        ylim = [10, 5000],
        yscale = :log10,
        ylabel = "POPS (cm⁻³)",
    )

    N_pops_smps = map(i -> sum(dN_dlog_Dp_DMA_count_0[kk, i] .* dlnDp[kk]), 1:length(tf))

    p1 = plot!(
        tf,
        N_pops_smps;
        color = :darkred,
        label = :none
    )
 
    p2 = plot(
        tf,
        N_count0;
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        yscale = :log10,
        ylim = [100, 50000],
        ylabel = "SMPS (cm⁻³)",
    )

    ii = (ccn_supersaturation .== 1.0) .& (denuder_flag .== false)

    p3 = plot(
        tf[ii],
        N_ccn0[ii];
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        color = :black,
        label = :none,
        yscale = :log10,
        ylim = [50, 5000],
        ylabel = "CCN (cm⁻³)",
        xlim = (ss, ee),
    )

    ii = (ccn_supersaturation .== 0.6) .& (denuder_flag .== false)

    p3 = plot!(tf[ii], N_ccn0[ii]; color = :darkred, label = :none)

    ii = (ccn_supersaturation .== 0.2) .& (denuder_flag .== false)

    p3 = plot!(tf[ii], N_ccn0[ii]; color = :darkgoldenrod, label = :none)

    NN = mapfoldl(
        i ->
            reverse(dN_dlog_Dp_DMA_count_2[:, i]) ./ maximum(dN_dlog_Dp_DMA_count_2[:, i]),
        hcat,
        1:length(tf),
    )

    p4 = heatmap(
        tf,
        reverse(diameter_mobility_midpoints[1:86]),
        NN[15:100, :];
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        ylabel = "D (nm)",
        ylim = (10, 400),
        yscale = :log10,
        cbar = :none,
        c = :jet,
    )

    NP = mapfoldl(
        i -> dNdlog_Dp_POPS[:, i] ./ maximum(dNdlog_Dp_POPS[:, i]),
        hcat,
        1:length(tf),
    )

    p5 = heatmap(
        tf,
        diameter_pops_midpoints,
        log10.(NP);
        xticks = (tm, ticks),
        minorgrid = true,
        minorticks = 6,
        ylabel = "D (nm)",
        ylim = (100, 4000),
        yscale = :log10,
        cbar = :none,
        c = :jet,
    )

    NC = mapfoldl(
        i ->
            droplet_size_distribution_ccn[:, i] ./
            maximum(droplet_size_distribution_ccn[:, i]),
        hcat,
        1:length(tf),
    )
    ND = map(
        i -> sum(
            collect(1:20) .* droplet_size_distribution_ccn[:, i] ./
            sum(droplet_size_distribution_ccn[:, i]),
        ),
        1:length(tf),
    )

    ii = (ccn_supersaturation .== 0.6) .& (denuder_flag .== false)

    p6 = scatter(
        tf[ii],
        ND[ii];
        xticks = (tm, ticks),
        xlim = (ss, ee),
        label = :none,
        ylabel = "DSD bin (-)",
        ylim = [7, 10],
        color = :black,
    )

    jj = (ccn_supersaturation .== 0.6) .& (denuder_flag .== true)

    p6 = scatter!(tf[jj], ND[jj]; color = :darkred, label = :none)

    mft = map(1:length(tf)) do i
        return fitme(
            ccn_supersaturation,
            diameter_mobility_midpoints,
            response_function_ccn,
            response_function_cpc_count,
            i,
        )
    end

    getme(x) =
        try
            #display(x[1])
            x[2][1]
            #sleep(0.1)
        catch
            plot([0, 0])
        end

    mp = map(getme, mft)

    ii = (denuder_flag .== false) .& (ccn_supersaturation .>= 0.4)

    p7 = scatter(
        tf[ii],
        mp[ii];
        xticks = (tm, ticks),
        xlim = (ss, ee),
        label = :none,
        ylabel = "κ (-)",
        color = :black,
        ylim = [0, 0.6],
    )

    jj = (denuder_flag .== true) .& (ccn_supersaturation .>= 0.4)

    p7 = scatter!(tf[jj], mp[jj]; label = :none, color = :darkred)

    ii = (denuder_flag .== false) .& (ccn_supersaturation .< 0.4)

    p8 = scatter(
        tf[ii],
        mp[ii];
        xticks = (tm, ticks),
        xlim = (ss, ee),
        label = :none,
        ylabel = "κ (-)",
        color = :black,
        ylim = [0, 1.2],
    )

    jj = (denuder_flag .== true) .& (ccn_supersaturation .< 0.4)

    p8 = scatter!(tf[jj], mp[jj]; label = :none, color = :darkred)

    pts = plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        p7,
        p8;
        xlim = (ss, ee),
        layout = grid(8, 1),
        size = (900, 1000),
        left_margin = 60px,
        right_margin = 20px,
    )

    savefig(pts, "src/assets/quicklooks/" * dstr * "_ql.png")

    return :SUCCESS
end

function my_quicklooks(date)
    try
        quicklook(date)
        
        open("src/quicklook.md", "a") do f
            d0 = Dates.format(date, "yyyymmdd")
            d1 = Dates.format(date, "mm/dd")
            d2 = Dates.format(date+Day(7), "mm/dd")
            write(f, "## $d1 - $d2\n")
            write(f, "\n")
            write(f, "![](assets/quicklooks/$(d0)_ql.png)\n\n")
        end

        open("src/housekeeping.md", "a") do f
            d0 = Dates.format(date, "yyyymmdd")
            d1 = Dates.format(date, "mm/dd")
            d2 = Dates.format(date+Day(7), "mm/dd")
            write(f, "## $d1 - $d2\n")
            write(f, "\n")
            write(f, "![](assets/hk/$(d0)_hk.png)\n\n")
        end

        :SUCCESS
    catch
        :FAILED
    end
end


