using NetCDF
using Dates

function read_netCDF(dp, ss)
    sss = ss:Day(1):ss+Day(7)
    dstr = Dates.format.(sss, "yyyymmdd")
 

    dfiles = map(x -> "epc10denudedccnS2.b1." * x * ".000000.cdf", dstr)
    avfiles = readdir(dp)
    files = intersect(dfiles, avfiles)
    
    readall(var, files; fun = vcat) = mapfoldl(f -> ncread(dp * f, var), fun, files)
    readfirst(var, files; fun = vcat) = mapfoldl(f -> ncread(dp * f, var), fun, [files[1]])

    ta = readall("base_time", files)
    total_N_conc_pops = readall("total_N_conc_pops", files)
    sample_flow_pops = readall("sample_flow_pops", files)
    sample_flow_dma = readall("sample_flow_dma", files)
    dNdlog_Dp_POPS = readall("dN_dlog_Dp_POPS", files; fun = hcat)
    diameter_pops_midpoints = readfirst("diameter_pops_midpoints", files)
    diameter_pops_bounds = readfirst("diameter_pops_bounds", files)
    sample_temperature = readall("sample_temperature", files)
    sample_relative_humidity = readall("sample_relative_humidity", files)
    inversion_matrix = readall("inversion_matrix", files)
    diameter_mobility_midpoints = readfirst("diameter_mobility_midpoints", files)
    diameter_mobility_bounds = readfirst("diameter_mobility_bounds", files)
    DMA_inner_radius = readfirst("DMA_inner_radius", files)
    DMA_outer_radius = readfirst("DMA_outer_radius", files)
    DMA_effective_length = readfirst("DMA_effective_length", files)
    hv_polarity = readfirst("hv_polarity", files)
    sheath_flow_dma = readall("sheath_flow_dma", files)
    sample_flow_dma = readall("sample_flow_dma", files)
    qc_cpc_flag = readall("qc_cpc_flag", files)
    ccn_gradient = readall("ccn_gradient", files)
    ccn_supersaturation = readall("ccn_supersaturation", files)
    denuder_flag = readall("denuder_flag", files)

    sample_flow_ccn = readall("sample_flow_ccn", files)
    sheath_flow_ccn = readall("sheath_flow_ccn", files)

    lat = readall("lat", files)
    lon = readall("lon", files)
    alt = readall("alt", files)

    response_function_cpc_serial = readall("response_function_cpc_serial", files; fun = hcat)
    response_function_cpc_count = readall("response_function_cpc_count", files; fun = hcat)
    response_function_ccn = readall("response_function_ccn", files; fun = hcat)
    dN_dlog_Dp_DMA_serial_0 = readall("dN_dlog_Dp_DMA_serial_0", files; fun = hcat)
    dN_dlog_Dp_DMA_serial_2 = readall("dN_dlog_Dp_DMA_serial_2", files; fun = hcat)
    dN_dlog_Dp_DMA_count_0 = readall("dN_dlog_Dp_DMA_count_0", files; fun = hcat)
    dN_dlog_Dp_DMA_count_2 = readall("dN_dlog_Dp_DMA_count_2", files; fun = hcat)
    dN_dlog_Dp_CCN_0 = readall("dN_dlog_Dp_CCN_0", files; fun = hcat)
    dN_dlog_Dp_CCN_2 = readall("dN_dlog_Dp_CCN_2", files; fun = hcat)
    droplet_size_distribution_ccn = readall("droplet_size_distribution_ccn", files; fun = hcat)

    tb = mapfoldl(i -> ncread(dp * files[i], "time_offset") .+ ta[i], vcat, 1:length(files))
    tf = map(unix2datetime, tb)

    return (
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
        sheath_flow_ccn
    )
end
