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
using Chain
using Plots

const conf = YAML.load_file("config.yaml")

const Vhi = Signal(150.0)
const Vlow = Signal(6000.0)
const tscan = Signal(120)
const thold = Signal(10)
const tflush = Signal(10)

include("dma_control.jl")               # High voltage power supply
include("labjack_io.jl")                # Labjack channels I/O
include("smps_signals.jl")              # Labjack channels I/O

const portPOPS = PrintedOpticalParticleSpectrometer.config(conf["serial"]["POPS"])
const portDMT = DropletMeasurementTechnologiesCCN.config(conf["serial"]["CCN"])
const portCPC = CondensationParticleCounters.config(Symbol(conf["CPC"]["model"]), conf["serial"]["CPC"])
const HANDLE = openUSBConnection(conf["LJ"]["ID"])
const U3HANDLE = u3.U3()
const caliInfo = getCalibrationInformation(HANDLE)
const caliInfoTdac = getTdacCalibrationInformation(HANDLE, conf["LJ"]["FIO"]["Tick_DAC"])
const Λ = get_DMA_config(conf["DMA"]["Qsh"], conf["DMA"]["Qsa"], conf["DMA"]["T"], conf["DMA"]["p"], Symbol(conf["DMA"]["model"]))

const calvolt = get_cal("voltage_calibration.csv")

const classifierV = Signal(200.0)
const dmaState = Signal(:SMPS)
const oneHz = every(1.0)
const smps_start_time = Signal(datetime2unix(now(UTC)))
const smps_elapsed_time = map(t -> Int(round(t - smps_start_time.value, digits=0)), oneHz)
const smps_scan_state, reset, V, Dp = smps_signals()
const signalV = map(calibrateVoltage, V)
const labjack_signals = map(labjackReadWrite, signalV)

function get_current_record()
    AIN, Tk, rawcount, count = labjack_signals.value
    RH = AIN[conf["LJ"]["AIN"]["RH"]+1] * 100.0
    T = AIN[conf["LJ"]["AIN"]["T"]+1] * 100.0 - 40.0
    readV = AIN[conf["LJ"]["AIN"]["V"]+1] |> (x -> (x * 1000.0))
    readI = AIN[conf["LJ"]["AIN"]["I"]+1] |> (x -> -x * 0.167 * 1000.0)
    @sprintf("LABJCACK,%i,%.3f,%s,%.3f,%.3f,%.3f,%.3f,%.3f", smps_elapsed_time.value, V.value, smps_scan_state.value, readV, readI, RH, T, count[1]./16.666666)
end

function start_acquisition_loops()
    @async PrintedOpticalParticleSpectrometer.stream(portPOPS, "pops.txt")
    @async DropletMeasurementTechnologiesCCN.stream(portDMT, "data.txt")
    @async CondensationParticleCounters.stream(portCPC, Symbol(conf["CPC"]["model"]), "cpc.txt")
end

function packet()
    cpc = CondensationParticleCounters.get_current_record()
    pops = PrintedOpticalParticleSpectrometer.get_current_record()
    ccn = DropletMeasurementTechnologiesCCN.get_current_record()
    lj = get_current_record()
    tc = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")

    return mapfoldl(x -> string(x) * ";", *, [tc, cpc, pops, ccn, lj])[1:end-1] * '\n'
end

start_acquisition_loops()

#U3HANDLE.configIO(EnableCounter0 = false, EnableCounter1 = false, NumberOfTimersEnabled = 0, FIOAnalog = 0)
#U3HANDLE.getFeedback(u3.BitStateWrite(10, 1))

file = "all.txt"
run(`touch $(file)`)

function acquire()
    x = packet()
    filter(x -> x != '\r', x)
    open(file, "a") do io
        write(io, x)
    end
end

asdf= map(x -> acquire(),oneHz)

const dataBufferdt = CircularBuffer{DateTime}(600)
const dataBufferCPCs = CircularBuffer{Float64}(600)
const dataBufferCPCc = CircularBuffer{Float64}(600)
const dataBufferCCN = CircularBuffer{Float64}(600)
const dataBufferVr = CircularBuffer{Float64}(600)
const dataBufferIr = CircularBuffer{Float64}(600)
packet()
function graphit()
    x = packet()
    a = split(x, ";")
    t = DateTime(a[1])
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
    del = @chain (dataBufferdt .- dataBufferdt[1]) Dates.value.(_) _./1000 
    p1 = plot(del, dataBufferCPCc, legend = false, color = :black)
    p1 = plot!(del, dataBufferCCN, legend = false, color = :darkred)
    p2 = plot(del, dataBufferVr, legend = false, color = :black, yscale = :log10, ylim = (100, 10000))
    #p3 = plot(del, dataBufferIr, legend = false, color = :black)

    pa = plot(p1, p2,  layout = grid(2,1))
    push!(p, pa)
end

p = Signal(plot([now()], [1.0]))

asd = map(_ -> graphit(), oneHz)
y = map(display, p)

function psd(ii, jj)
    Dd = map(v -> vtod(Λ, v), dataBufferVr)
    N = map(x -> x, dataBufferCPCc)
    Nccn = map(x -> x, dataBufferCCN)
    p = plot(Dd[1:end-ii], N[ii+1:end], xscale = :log10, xlim = [40, 1000], ylim = [0,100], xticks = [100, 1000], minorgrid = true)
    p = plot!(Dd[1:end-jj], Nccn[jj+1:end])
    mDe = exp10.(range(log10(40), stop = log10(400), length = 101))
    mDd = sqrt.(mDe[2:end] .* mDe[1:end-1])

    #itp1 = interpolate(Dd[1:end-ii], N[ii+1:end])
    display(p)
end

psd(4, 24)

display(p.value)