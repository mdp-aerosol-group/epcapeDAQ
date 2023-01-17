using DropletMeasurementTechnologiesCCN
using CondensationParticleCounters
using PrintedOpticalParticleSpectrometer
using DifferentialMobilityAnalyzers
using LabjackU6Library

using PyCall
using CSV
using YAML
using DataFrames
using Reactive
using Dates
using Printf
using Interpolations
using DataStructures
using CircularList
using LibSerialPort
using Chain
using Plots
unicodeplots()

const bp = "/home/aerosol/Data/"
const conf = YAML.load_file("config.yaml")

struct State
    dT::Float64     # CCN temperature gradient
    denuded::Bool   # Denuded/undenuded
end

state = circularlist(State(4.0, false))
insert!(state, State(4.0, true))
insert!(state, State(6.0, false))
insert!(state, State(6.0, true))
insert!(state, State(10.0, false))
insert!(state, State(10.0, true))
insert!(state, State(14.0, false))
insert!(state, State(14.0, true))

const Vhi = Signal(150.0)
const Vlow = Signal(6000.0)
const tscan = Signal(240)
const thold = Signal(100)
const tflush = Signal(10)

const u3 = pyimport("u3")
const U3HANDLE = u3.U3()

include("dma_control.jl")               # High voltage power supply
include("labjack_io.jl")                # Labjack channels I/O
include("smps_signals.jl")              # Labjack channels I/O
include("valworx.jl")                   # Labjack channels I/O

const portPOPS = PrintedOpticalParticleSpectrometer.config(conf["serial"]["POPS"])
const portDMT = DropletMeasurementTechnologiesCCN.config(conf["serial"]["CCN"])
const portCPC =
    CondensationParticleCounters.config(Symbol(conf["CPC"]["model"]), conf["serial"]["CPC"])
const HANDLE = openUSBConnection(conf["LJ"]["ID"])
const caliInfo = getCalibrationInformation(HANDLE)
const caliInfoTdac = getTdacCalibrationInformation(HANDLE, conf["LJ"]["FIO"]["Tick_DAC"])
const Λ = get_DMA_config(
    conf["DMA"]["Qsh"],
    conf["DMA"]["Qsa"],
    conf["DMA"]["T"],
    conf["DMA"]["p"],
    Symbol(conf["DMA"]["model"]),
)

const calvolt = get_cal("voltage_calibration.csv")

const classifierV = Signal(200.0)
const dmaState = Signal(:SMPS)
const oneHz = every(1.0)
const smps_start_time = Signal(datetime2unix(now(UTC)))
const smps_elapsed_time = map(t -> Int(round(t - smps_start_time.value; digits = 0)), oneHz)
const smps_scan_state, reset, V, Dp = smps_signals()
const signalV = map(calibrateVoltage, V)
const labjack_signals = map(labjackReadWrite, signalV)

function updatestate!()
    CircularList.forward!(state)
    DropletMeasurementTechnologiesCCN.set_dT(portDMT, CircularList.current(state).data.dT)
    return valve(CircularList.current(state).data.denuded)
end

stateLoop = map(_ -> updatestate!(), reset)

function get_current_record()
    AIN, Tk, rawcount, count = labjack_signals.value
    RH = AIN[conf["LJ"]["AIN"]["RH"]+1] * 100.0
    T = AIN[conf["LJ"]["AIN"]["T"]+1] * 100.0 - 40.0
    readV = AIN[conf["LJ"]["AIN"]["V"]+1] |> (x -> (x * 1000.0))
    readI = AIN[conf["LJ"]["AIN"]["I"]+1] |> (x -> -x * 0.167 * 1000.0)
    @sprintf(
        "LABJCACK,%i,%.3f,%s,%.3f,%.3f,%.3f,%.3f,%.3f",
        smps_elapsed_time.value,
        V.value,
        smps_scan_state.value,
        readV,
        readI,
        RH,
        T,
        count[1] ./ 16.666666
    )
end

function start_acquisition_loops()
    @async PrintedOpticalParticleSpectrometer.stream(portPOPS, bp * "mtsncsupopsrs232/pops")
    @async DropletMeasurementTechnologiesCCN.stream(portDMT, bp * "mtsncsuccn100/ccn")
    @async CondensationParticleCounters.stream(
        portCPC,
        Symbol(conf["CPC"]["model"]),
        bp * "mtsncsucpc3772/cpc",
    )
end

function packet()
    cpc = CondensationParticleCounters.get_current_record()
    pops = PrintedOpticalParticleSpectrometer.get_current_record()
    ccn = DropletMeasurementTechnologiesCCN.get_current_record()
    lj = get_current_record()
    tc = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS") * "," * string(CircularList.current(state).data.denuded)

    return mapfoldl(x -> string(x) * ";", *, [tc, cpc, pops, ccn, lj])[1:end-1] * '\n'
end

start_acquisition_loops()

function acquire()
    x = packet()
    filter(x -> x != '\r', x)
    tc = Dates.format(now(), "yyyymmdd")
    open(bp * "mtsncsudenudedccn/rack" * "_" * tc * ".txt", "a") do io
        return write(io, x)
    end
end

daqLoop = map(x -> acquire(), oneHz)

const dataBufferdt = CircularBuffer{DateTime}(600)
const dataBufferCPCs = CircularBuffer{Float64}(600)
const dataBufferCPCc = CircularBuffer{Float64}(600)
const dataBufferCCN = CircularBuffer{Float64}(600)
const dataBufferVr = CircularBuffer{Float64}(600)
const dataBufferIr = CircularBuffer{Float64}(600)

function graphit()
    x = packet()
    a = split(x, ";") 
    b = split(a[1], ",")
    t = DateTime(b[1])
    cs = @chain split(a[2], ",") getindex(_, 3) parse(Float64, _)
    cc = @chain split(a[5], ",") getindex(_, 9) parse(Float64, _)
    Vr = @chain split(a[5], ",") getindex(_, 5) parse(Float64, _)
    Ir = @chain split(a[5], ",") getindex(_, 6) parse(Float64, _)
    ccn = @chain split(a[4], ",") getindex(_, 21) parse(Float64, _)
    push!(dataBufferdt, t)
    push!(dataBufferCPCs, cs)
    push!(dataBufferCPCc, cc)
    push!(dataBufferVr, Vr)
    push!(dataBufferIr, Ir)
    push!(dataBufferCCN, ccn)
    del = @chain (dataBufferdt .- dataBufferdt[1]) Dates.value.(_) _ ./ 1000
    p1 = plot(del, dataBufferCPCc; legend = false, color = :black)
    p1 = plot!(del, dataBufferCCN; legend = false, color = :darkred)
    p2 = plot(
        del,
        dataBufferVr;
        legend = false,
        color = :black,
        yscale = :log10,
        ylim = (100, 10000),
    )

    pa = plot(p1, p2; layout = grid(2, 1))
    if typeof(p0) == Plots.Plot{Plots.UnicodePlotsBackend}
        push!(p, p1)
    end
end

p0 = plot([now()], [1.0])
p = Signal(p0)

sleep(4)
graphLoop = map(_ -> graphit(), oneHz)

graphHz = every(60.0)
graphDisp = map(graphHz) do _
    display(packet())
    display(p.value)
end

# function psd(ii, jj)
#     Dd = map(v -> vtod(Λ, v), dataBufferVr)
#     N = map(x -> x, dataBufferCPCc)
#     Nccn = map(x -> x, dataBufferCCN)
#     p = plot(
#         Dd[1:end-ii],
#         N[ii+1:end];
#         xscale = :log10,
#         xlim = [40, 1000],
#         ylim = [0, 100],
#         xticks = [100, 1000],
#         minorgrid = true,
#     )
#     p = plot!(Dd[1:end-jj], Nccn[jj+1:end])
#     mDe = exp10.(range(log10(40); stop = log10(400), length = 101))
#     mDd = sqrt.(mDe[2:end] .* mDe[1:end-1])

#     #itp1 = interpolate(Dd[1:end-ii], N[ii+1:end])
#     return display(p)
# end

# psd(4, 24)

# display(p.value)