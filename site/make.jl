using Documenter
using DelimitedFiles
using Chain
using Plots
using Dates
using TimeSeries
using DataFrames
using Statistics
using MLStyle
using DifferentialMobilityAnalyzers

include("processingroutines.jl")
include("cdf_header.jl")
include("quicklooks.jl")

function make_netcdf(ss)
    cdfname = create_file(ss)
    t, denuded, tcpc, tpops, tccn, tlj = readraw(ss)
    Λ = get_DMA_config(5.0, 1.5, 20.0, 101315.0, :TSILONG)
    δ = setupDMA(Λ, vtoz(Λ, 7000), vtoz(Λ, 30), 100)
    #δ = setupSMPS(Λ, 30, 7000, 240, 2.4)

    time2cdf(cdfname, t, ss)
    pops2cdf(cdfname, t, tpops, ss)
    state2cdf(cdfname, t, tlj, ss)
    dmasetup2cdf(cdfname, Λ, δ)
    ccnflow2cdf(cdfname, t, tccn, ss)
    distribution2cdf(cdfname, t, tcpc, tccn, tlj, denuded, Λ, δ, ss)
    latlonalt2cdf(cdfname)

    return nothing
end

make_netcdf(today()+Day(1))
today() |> make_netcdf

dates = Date(2023, 02, 10):Day(7):Date(2023, 08, 31)

try
    rm("src/quicklook.md")
    touch("src/quicklook.md")
    rm("src/housekeeping.md")
    touch("src/housekeeping.md")
catch
end

map(my_quicklooks, dates)

makedocs(;
    sitename = "EPCAPE-CCN",
    authors = "Markus Petters",
    pages = Any[
        "Home"=>"index.md",
        "Quicklooks"=>"quicklook.md",
        "Housekeeping"=>"housekeeping.md",
    ],
)
