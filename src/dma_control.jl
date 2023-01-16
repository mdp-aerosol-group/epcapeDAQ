function get_cal(file)
    df = CSV.read(file, DataFrame)
    itp = interpolate((df[!, :readV],), df[!, :setV], Gridded(Linear()))
    extp = extrapolate(itp, Flat())

    return extp
end
vtod(Λ, v) = @chain vtoz(Λ, v) ztod(Λ, 1, _)


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
    qsa = Qsa * qsh
    t = T + 273.15
    leff = 0.0
    polarity = :+

    Λ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, polarity, 6, form)

    return Λ
end

calibrateVoltage(v) = getVdac(calvolt(v), :+, true)

function getVdac(setV::Float64, polarity::Symbol, powerSwitch::Bool)
    (setV > 0.0) || (setV = 0.0)
    (setV < 10000.0) || (setV = 10000.0)

    if polarity == :-
        # Negative power supply +0.36V = -10kV, 5V = 0kV
        m = 10000.0 / (0.36 - 5.03)
        b = 10000.0 - m * 0.36
        setVdac = (setV - b) / m
        if setVdac < 0.36
            setVdac = 0.36
        elseif setVdac > 5.1
            setVdac = 5.1
        end
        if powerSwitch == false
            setVdac = 5.0
        end
    elseif polarity == :+
        # Positive power supply +0V = 0kV, 4.64V = 0kV
        m = 10000.0 / (4.64 - 0)
        b = 0
        setVdac = (setV - b) / m
        if setVdac < 0.0
            setVdac = 0.0
        elseif setVdac > 4.64
            setVdac = 4.64
        end
        if powerSwitch == false
            setVdac = 0.0
        end
    end
    return setVdac
end