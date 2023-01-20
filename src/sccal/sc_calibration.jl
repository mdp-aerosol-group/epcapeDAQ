using SpecialFunctions
using LsqFit
using MLStyle
using Interpolations
using NetCDF
using Dates
using Plots

function threshold(N, Dp, c::Float64, n1::Float64, n2::Float64)
    M = deepcopy(N)
    M[(N.<=c).&(Dp.>150)] .= n2
    M[(N.<=c).&(Dp.<150)] .= n1
    return M
end

Taf(Dp, μ, σ) = @. 0.5 * (1.0 + erf((Dp - μ) ./ (sqrt(2.0 .* σ))))
model(x, p) = p[3] * Taf(x, p[1], p[2])
linfit(x, y) = [ones(length(x), 1) x] \ y

ig(x) = @match x begin
    4.0 => 100.0
    6.0 => 58.0
    8.0 => 47.0
    10.0 => 40.0
    12.0 => 35.0
    14.0 => 30.0
    16.0 => 27.0
    _ => 50.0
end

# Calibration at NCSU with dried ammonium sulfate.
dp = ""
files = ["epc10denudedccnA1.b1.20230119.000000.cdf"]

readall(var, files) = mapfoldl(f -> ncread(dp * f, var), vcat, files)

ta = readall("base_time", files)
total_N_conc_pops = readall("total_N_conc_pops", files)
sample_flow_pops = readall("sample_flow_pops", files)
sample_flow_dma = readall("sample_flow_dma", files)
dNdlog_Dp_POPS = readall("dN_dlog_Dp_POPS", files)
diameter_pops_midpoints = readall("diameter_pops_midpoints", files)
diameter_pops_bounds = readall("diameter_pops_bounds", files)
sample_temperature = readall("sample_temperature", files)
sample_relative_humidity = readall("sample_relative_humidity", files)
inversion_matrix = readall("inversion_matrix", files)
diameter_mobility_midpoints = readall("diameter_mobility_midpoints", files)
diameter_mobility_bounds = readall("diameter_mobility_bounds", files)
DMA_inner_radius = readall("DMA_inner_radius", files)
DMA_outer_radius = readall("DMA_outer_radius", files)
DMA_effective_length = readall("DMA_effective_length", files)
hv_polarity = readall("hv_polarity", files)
sheath_flow_dma = readall("sheath_flow_dma", files)
sample_flow_dma = readall("sample_flow_dma", files)
qc_cpc_flag = readall("qc_cpc_flag", files)
ccn_gradient = readall("ccn_gradient", files)
denuder_flag = readall("denuder_flag", files)

response_function_cpc_serial = readall("response_function_cpc_serial", files)
response_function_cpc_count = readall("response_function_cpc_count", files)
response_function_ccn = readall("response_function_ccn", files)
dN_dlog_Dp_DMA_serial_0 = readall("dN_dlog_Dp_DMA_serial_0", files)
dN_dlog_Dp_DMA_serial_2 = readall("dN_dlog_Dp_DMA_serial_2", files)
dN_dlog_Dp_DMA_count_0 = readall("dN_dlog_Dp_DMA_count_0", files)
dN_dlog_Dp_DMA_count_2 = readall("dN_dlog_Dp_DMA_count_2", files)
dN_dlog_Dp_CCN_0 = readall("dN_dlog_Dp_CCN_0", files)
dN_dlog_Dp_CCN_2 = readall("dN_dlog_Dp_CCN_2", files)
droplet_size_distribution_ccn = readall("droplet_size_distribution_ccn", files)

tb = mapfoldl(i -> ncread(dp * files[i], "time_offset") .+ ta[i], vcat, 1:length(files))
tf = map(unix2datetime, tb)

include("rh_functions.jl")
include("")

ii = tf .> DateTime(2023, 1, 19, 15, 0, 0)
rCCN = response_function_ccn[:, ii]
rCN = response_function_cpc_serial[:, ii]
nCCN = dN_dlog_Dp_CCN_0[:, ii]
nCN = dN_dlog_Dp_DMA_count_0[:, ii]
n2CCN = dN_dlog_Dp_CCN_2[:, ii]
n2CN = dN_dlog_Dp_DMA_count_2[:, ii]

isdenuded = denuder_flag[ii]
Dp = diameter_mobility_midpoints
dT = ccn_gradient[ii]
i = 8

blist = [
    [4, 10, 11, 14, 15, 18, 19, 20, 40, 43, 45, 46, 47, 48, 67, 71, 74, 75]
    collect(89:107)
]

foo = map(1:length(dT)) do i
    if i ∈ blist
        return nothing
    end

    p1a = plot(
        Dp,
        nCN[:, i];
        xscale = :log10,
        label = "rCN",
        xlim = (10, 400),
        minorgrid = true,
    )
    p1a = plot!(Dp, nCCN[:, i]; label = "rCCN: $(isdenuded[i])")

    p2a = scatter(
        Dp,
        nCCN[:, i] ./ nCN[:, i];
        ylim = (0, 1.4),
        xscale = :log10,
        xlim = (10, 400),
        minorgrid = true,
        label = "rCCN/rCN $(dT[i])",
        legend = :topleft,
    )

    tCCN = threshold(nCCN[:, i], Dp, 0.1, 0.1, 0.1)
    tCN = threshold(nCN[:, i], Dp, 0.1, 0.0, 0.1)
    af = tCCN ./ tCN
    jj = (Dp .< 180.0) .& isfinite.(af)

    initial = [ig(dT[i]), ig(dT[i]) ./ 10, 0.9]
    Ax = try
        fit = curve_fit(model, Dp[jj], af[jj], initial)
        if fit.converged == true
            fit.param
        else
            return nothing
        end
    catch
        return nothing
    end

    plot!(Dp, model(Dp, Ax); label = "$(round(Ax[1]))")

    pp1a = plot(p1a, p2a; layout = grid(2, 1), title = "$i")

    println(dT[i], " ", Ax[1], " ", isdenuded[i])
    display(pp1a)
    #sleep(1)
    if Ax[1] .<= 0
        return nothing
    end
    return dT[i], Ax[1], isdenuded[i]
end

extract = filter(x -> ~isnothing(x), foo)

mydT = map(x -> x[1], extract)
myD50 = map(x -> x[2], extract)
myDen = map(x -> x[3], extract)

ii = myDen .== 0
jj = myDen .== 1

mdT = 4:1:16
thesc = get_sc.(myD50 .* 1e-9)
mfit = linfit(mydT, thesc)

p = scatter(
    mydT[jj],
    thesc[jj];
    legend = :topleft,
    alpha = 0.4,
    ylim = [0, 1.4],
    label = "denuded",
    xlabel = "dT (K)",
    ylabel = "Superaturation (%)",
)

println(mfit)
m = round(mfit[2]; digits = 4)
b = round(mfit[1]; digits = 5)

p = scatter!(mydT[ii], thesc[ii]; alpha = 0.4, label = "undenuded")
p = plot!(mdT, m .* mdT .+ b; color = :black, label = "calibration")
x0, x1 = xlims(p)
y0, y1 = ylims(p)
dx = x1 - x0
dy = y1 - y0
annotate!(x0 + 0.04 * dx, y1 - 0.23 * dy, text("sc = $m dT + $b", :black, :left, 10))

savefig("sc_calibration.pdf")