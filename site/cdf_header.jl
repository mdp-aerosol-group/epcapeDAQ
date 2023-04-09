using NetCDF
using Dates
using Chain

function create_file(ss)
    globalatts = Dict(
        "site_id" => "epc",
        "platform_id" => "denudedccn",
        "facility_id" => "S2",
        "data_level" => "b1",
        "location_description" => "Eastern Pacific Cloud Aerosol Precipitation Experiment (EPCAPE), Mount Soledad",
        "datastream" => "epc10denudedccnS2.b1",
        "sampling_interval" => "10 min",
        "hold" => "100 s",
        "up_scan" => "240 s",
        "up_hold" => "10 s",
        "down_hold" => "240 s",
        "flush" => "10 s",
        "DMA_model" => "TSI 3071",
        "neutralizer_model" => "ADI 100, 210Po",
        "detector_model_cn" => "TSI 3772",
        "detector_model_ccn" => "DMT 100",
        "CCN duty cycle" => "ss = 0.2, 0.4, 0.6, 0.8, 1.0%",
        "detector_sample_flow_cn" => "1 L min-1",
        "instrument_output_units_cn" => "dN/dlogDp (natural log)",
        "detector_sample_flow_ccn" => "0.3 L min-1",
        "instrument_output_units_ccn" => "dN/dlogDp (natural log)",
        "auxiliary_detector_model_opc" => "POPS",
        "auxiliary_detector_sample_flow_opc" => "0.3 L min-1",
        "instrument_output_units_opc" => "dN/dlogDp (natural log)",

        "history" => "created by NCSU, contact Markus Petters: mdpetter@ncsu.edu",
    )

    datestr = Dates.format(DateTime(ss), "yyyymmdd.HHMMSS")
    filename = "/home/mdpetter/epcape/epc10denudedccnS2.b1/epc10denudedccnS2.b1." * datestr * ".cdf"

    dtav = Minute(5)
    ts = DateTime(ss):dtav:(DateTime(ss)+Day(1)-Second(1))

    # Dimensions
    timedim = @chain datetime2unix.(ts) collect(_) Int32.(_)
    mobdim = zeros(100)
    opcdim = zeros(16)
    bound = zeros(Int32, 2)
    ccnopcdim = zeros(20)
    dutydim = zeros(600)

    # Dimension attributes
    mobatts =
        Dict("longname" => "Midpoint of geometric mean mobility diameter", "units" => "nm")
    opcatts =
        Dict("longname" => "Midpoint of PSL equivalent optical diameter", "units" => "nm")
    timatts = Dict(
        "longname" => "Base time in Epoch",
        "string" => Dates.format(DateTime(ss), "yyyy-mm-dd HH:MM:SS"),
        "units" => "seconds since 1970-1-1 0:00:00 0:00",
        "ancillary_variables" => "time_offset",
    )

    isfile(filename) && rm(filename)

    # Create Variables

    varatts = Dict(
        "longname" => "Base time in Epoch",
        "string" => Dates.format(DateTime(ss), "yyyy-mm-dd HH:MM:SS"),
        "units" => "seconds since 1970-1-1 0:00:00 0:00",
        "ancillary_variables" => "time_offset",
    )

    nccreate(filename, "base_time"; atts = varatts, t = NC_INT)

    varatts = Dict(
        "longname" => "Time offset from base_time",
        "string" => Dates.format(DateTime(ss), "yyyy-mm-dd HH:MM:SS"),
        "units" => "seconds since 1970-1-1 0:00:00 0:00",
        "ancillary_variables" => "time_offset",
    )

    nccreate(filename, "time_offset", "time", timedim, timatts; atts = varatts, t = NC_INT)

    varatts = Dict(
        "longname" => "Number size distribution, electrical mobility diameter derived from serial signal.",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "0th order Tikhonov using an initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The initial guess is method is very robust. However, measurement noise enters into the initial guess, which leads to a noise inverted distribution. ",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_DMA_serial_0",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "DMA response function derived from the serial signal of the CPC. ",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "The response function is the uninverted concentration measured by the detector after the DMA.",
    )

    nccreate(
        filename,
        "response_function_cpc_serial",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "DMA response function derived from the pulse counter of the CPC. ",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "The response function is the uninverted concentration measured by the detector after the DMA.",
    )

    nccreate(
        filename,
        "response_function_cpc_count",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "DMA response function derived from the concentration of the CCN. ",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "The response function is the uninverted concentration measured by the detector after the DMA.",
    )

    nccreate(
        filename,
        "response_function_ccn",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "Number size distribution, electrical mobility diameter derived from serial signal.",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "0th order Tikhonov using an initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The initial guess is method is very robust. However, measurement noise enters into the initial guess, which leads to a noise inverted distribution. ",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_DMA_count_0",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "CCN size distribution",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "0th order Tikhonov using an initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The initial guess is method is very robust. However, measurement noise enters into the initial guess, which leads to a noise inverted distribution. ",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_CCN_0",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )


    varatts = Dict(
        "longname" => "Number size distribution, electrical mobility diameter derived from serial signal.",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "2nd order Tikhonov using no initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The 2nd order Tikhonov inversion yields a smoothed solution that removes measurement noise. Compare with dN_dlog_Dp_DMA0 to ensure that the solution is not oversmoothed",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_DMA_serial_2",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "Number size distribution, electrical mobility diameter derived from pulse counter",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "2nd order Tikhonov using no initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The 2nd order Tikhonov inversion yields a smoothed solution that removes measurement noise. Compare with dN_dlog_Dp_DMA0 to ensure that the solution is not oversmoothed",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_DMA_count_2",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )


    varatts = Dict(
        "longname" => "CCN size distribution",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
        "inversion_method" => "2nd order Tikhonov using no initial guess",
        "inversion_software" => "DifferentialMobilityAnalyzers.jl v2.5.6 (https://doi.org/10.5194/amt-14-7909-2021)",
        "inversion_comment" => "The 2nd order Tikhonov inversion yields a smoothed solution that removes measurement noise. Compare with dN_dlog_Dp_DMA0 to ensure that the solution is not oversmoothed",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_CCN_2",
        "diameter_mobility",
        mobdim,
        mobatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    # OPC -----------------------------------------------------------------------------------------------------

    varatts = Dict(
        "longname" => "Total number concentration from integrated size distribution, POPS",
        "units" => "cm^3",
        "missing_value" => "NaN",
    )

    nccreate(filename, "total_N_conc_pops", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Sample flow rate, POPS",
        "units" => "cm^3",
        "missing_value" => "NaN",
    )

    nccreate(filename, "sample_flow_pops", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Number size distribution, printed optical spectrometer",
        "units" => "1/cm^3",
        "missing_value" => "NaN",
        "comment" => "dN_dlogDp is the aerosol number size distribution where the number of particles per bin (dN) have been divided by the bin-width in natural log space (dlogDp).  This simplifies comparison of size distributions from instruments with different bin spacing.",
    )

    nccreate(
        filename,
        "dN_dlog_Dp_POPS",
        "diameter_opc",
        opcdim,
        opcatts,
        "time",
        timedim,
        timatts;
        atts = varatts,
    )

    varatts = Dict("longname" => "POPS diameter bin boundaries", "units" => "nm", "comment" => "Based on calibration with size-selected ammonium sulfate n = 1.55 + 0i and a scaled Mie model")

    nccreate(
        filename,
        "diameter_pops_bounds",
        "bound",
        bound,
        "diameter_opc",
        opcdim,
        opcatts;
        atts = varatts,
    )

    varatts = Dict("longname" => "POPS diameter bin Midpoints", "units" => "nm")

    nccreate(
        filename,
        "diameter_pops_midpoints",
        "diameter_opc",
        opcdim,
        opcatts;
        atts = varatts,
    )

    # Temperature and RH sensor ------------------------------------------------------------------

    varatts = Dict(
        "longname" => "Relative humidty",
        "units" => "%",
        "missing_value" => "NaN",
        "comment" => "Rotronic HC2 temperature and humidity probe measuring at the inlet of the DMA system before the neutralizer.",
    )

    nccreate(filename, "sample_relative_humidity", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Sample temperature",
        "units" => "degC",
        "missing_value" => "NaN",
        "comment" => "Rotronic HC2 temperature and humidity probe measuring at the inlet of the DMA system before the neutralizer.",
    )

    nccreate(filename, "sample_temperature", "time", timedim, timatts; atts = varatts)


    varatts = Dict(
        "longname" => "Sample pressure",
        "units" => "kPa",
        "missing_value" => "NaN",
        "comment" => "Absolute pressure measured by the CPC at the exit of the DMA",
    )

    nccreate(filename, "sample_pressure", "time", timedim, timatts; atts = varatts)

    # Location ------------------------------------------------------------------

    varatts = Dict(
        "longname" => "Altitude above mean sea level",
        "units" => "m",
        "standard_name" => "altitude",
    )

    nccreate(filename, "alt"; atts = varatts)

    varatts = Dict(
        "longname" => "North latitude",
        "units" => "degree_N",
        "valid_min" => "-90.f",
        "valid_max" => "90.f",
        "standard_name" => "latitude",
    )

    nccreate(filename, "lat"; atts = varatts)

    varatts = Dict(
        "longname" => "East longitude",
        "units" => "degree_E",
        "valid_min" => "-180.f",
        "valid_max" => "180.f",
        "standard_name" => "longitude",
    )

    nccreate(filename, "lon"; atts = varatts)

    # DMA System ----------------------------------------------------------------------------------------

    varatts = Dict(
        "longname" => "DMA sample flow rate (CPC + CCN)",
        "units" => "m^3/s^1",
        "missing_value" => "NaN",
    )

    nccreate(filename, "sample_flow_dma"; atts = varatts)

    varatts = Dict(
        "longname" => "DMA sheath flow rate set point",
        "units" => "m^3/s^1",
        "missing_value" => "NaN",
    )

    nccreate(filename, "sheath_flow_dma"; atts = varatts)

    varatts = Dict(
        "longname" => "Time required for aerosol to flow from classifier exit to detector in CPC",
        "units" => "s",
        "missing_value" => "NaN",
    )

    nccreate(filename, "delay_time_cpc", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Time required for aerosol to flow from classifier exit to detector in CCN without denuder in line",
        "units" => "s",
        "missing_value" => "NaN",
    )

    nccreate(filename, "delay_time_ccn_undenuded", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Time required for aerosol to flow from classifier exit to detector in CCN with denuder in line",
        "units" => "s",
        "missing_value" => "NaN",
    )

    nccreate(filename, "delay_time_ccn_denuded", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Reference gas temperature",
        "units" => "K",
        "missing_value" => "NaN",
        "comment" => "Temperature used to compute the inversion matrix.",
    )

    nccreate(filename, "reference_gas_temperature"; atts = varatts)

    varatts = Dict(
        "longname" => "Reference gas pressure",
        "units" => "Pa",
        "missing_value" => "NaN",
        "comment" => "Pressure used to compute the inversion matrix.",
    )

    nccreate(filename, "reference_gas_pressure"; atts = varatts)

    varatts = Dict(
        "longname" => "HV polarity",
        "units" => "1",
        "missing_value" => "NaN",
        "flag_values" => "-1, 1",
        "flag_meanings" => "negative positive",
    )

    nccreate(filename, "hv_polarity"; atts = varatts)

    varatts =
        Dict("longname" => "Inner radius of DMA", "units" => "m", "missing_value" => "NaN")

    nccreate(filename, "DMA_inner_radius"; atts = varatts)

    varatts =
        Dict("longname" => "Outer radius of DMA", "units" => "m", "missing_value" => "NaN")

    nccreate(filename, "DMA_outer_radius"; atts = varatts)

    varatts = Dict("longname" => "Length of DMA", "units" => "m", "missing_value" => "NaN")

    nccreate(filename, "DMA_length"; atts = varatts)

    varatts = Dict(
        "longname" => "Effective length of DMA",
        "units" => "m",
        "missing_value" => "NaN",
        "comment" => "Parameterized diffusional losses inside the DMA entering the inversion matrix",
    )

    nccreate(filename, "DMA_effective_length"; atts = varatts)

    varatts = Dict("longname" => "Mobility diameter bin boundaries", "units" => "nm")

    nccreate(
        filename,
        "diameter_mobility_bounds",
        "bound",
        bound,
        "diameter_mobility",
        mobdim,
        mobatts;
        atts = varatts,
    )

    varatts = Dict("longname" => "DMA diameter bin Midpoints", "units" => "nm")

    nccreate(
        filename,
        "diameter_mobility_midpoints",
        "diameter_mobility",
        mobdim,
        mobatts;
        atts = varatts,
    )

    varatts = Dict(
        "longname" => "DMA Inversion matrix",
        "units" => "-",
        "missing value" => "NaN",
        "comment" => "Matrix used to convert response function to size distribution. The inversion accounts for multi-charge correction, diffusional broadening of the transfer function, and losses inside the DMA. See Petters (2018), https://doi.org/10.1080/02786826.2018.1530724 for details.",
    )

    nccreate(
        filename,
        "inversion_matrix",
        "diameter_mobility",
        mobdim,
        mobatts,
        "diameter_mobility",
        mobdim,
        mobatts;
        atts = varatts,
    )

    # ccn ----------------------------------------

    varatts = Dict(
        "longname" => "CCN sheath flow",
        "units" => "ml/min^1",
        "missing_value" => "NaN",
    )

    nccreate(filename, "sheath_flow_ccn", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "CCN sample flow",
        "units" => "ml/min^1",
        "missing_value" => "NaN",
    )

    nccreate(filename, "sample_flow_ccn", "time", timedim, timatts; atts = varatts)

    # QC ----------------------------------------

    varatts = Dict(
        "longname" => "QC flag for condensation particle counter",
        "units" => "16-bit integer in hexadecimal format",
        "missing value" => "NaN",
        "comment" => "RIE output from CPC. The parameter is in error if the bit is set. Bit 0x0001 => Saturator Temp, Bit 0x0002 => Condenser Temp, Bit 0x0004 => Optics Temp
, Bit 0x0008 => Inlet Flow Rate, Bit 0x0010 => Aerosol Flow Rate, Bit 0x0020 => Laser Power, Bit 0x0040 => Liquid Level, Bit 0x0080 => Concentration, Bit 0x0100 => Calibration Reminder. Higher bits are unused",
    )

    nccreate(filename, "qc_cpc_flag", "time", timedim, timatts; atts = varatts, t = NC_INT)


    # CCN parameters
    varatts = Dict(
        "longname" => "Nominal CCN thermal gradient",
        "units" => "K",
        "missing value" => "NaN",
        "comment" => "The streamwise thermal gradient sets the supersaturation in the CCN column",
    )

    nccreate(filename, "ccn_gradient", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Nominal CCN supersaturation",
        "units" => "%",
        "missing value" => "NaN",
        "calibration" => "Ammonium sulfate 20230119, ss = 0.0777*dT - 0.14895, dT in [K]",
        "related field" => "ccn_gradient",
        "comment" => "The supersaturation in the CCN column derived from thermal gradient and dried ammonium sulfate calibration",
    )

    nccreate(filename, "ccn_supersaturation", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Flag to indicate whether denuder was used prior to CCN measurement",
        "units" => "0 (no denuder) or 1 (denuder)",
        "missing value" => "NaN",
        "comment" => "The streamwise thermal gradient sets the supersaturation in the CCN column",
    )

    nccreate(filename, "denuder_flag", "time", timedim, timatts; atts = varatts, t = NC_INT)

    varatts = Dict(
        "longname" => "Flag to indicate whether the instrument was sampling from CVI",
        "units" => "0 to 1",
        "missing value" => "NaN",
        "comment" => "fraction of seconds if CVI valve was turned on, 1 = 100% on CVI",
    )

    nccreate(filename, "cvi_flag", "time", timedim, timatts; atts = varatts)

    varatts = Dict(
        "longname" => "Mean count distribution, DMT CCN optical spectrometer",
        "units" => "-",
        "missing_value" => "NaN",
        "comment" => "The DMT OPC measures the droplet size distribution at the exit of the column in counts/s. This field gives the average of this distribution during each scan.",
    )

    nccreate(
        filename,
        "droplet_size_distribution_ccn",
        "bin_ccn_opc",
        ccnopcdim,
        "time",
        timedim,
        timatts;
        atts = varatts,
    ) 
    ccnopcdim

    NetCDF.open(filename; mode = NC_WRITE) do f
        return NetCDF.putatt(f, "global", globalatts)
    end

    return filename
end

ss = Date(2022, 1, 1)

create_file(ss)