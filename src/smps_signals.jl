function smps_signals()
    holdTime, scanTime, flushTime, scanLength, startVoltage, endVoltage, c =
        scan_parameters()
    
    function trianglewave(x, T, l, u)
        modx = (mod(x, T) + T / 2)
        l + (u - l) * ifelse(modx < T, 2modx / T - 1, -2modx / T + 3)
    end

    # Set SMPS states
    function state(currentTime)
        holdTime, scanTime, flushTime, scanLength, startVoltage, endVoltage, c =
            scan_parameters()

        scanState = "DONE"
        (currentTime <= scanLength) && (scanState = "FLUSH")
        (currentTime < 2*scanTime + holdTime + flushTime) && (scanState = "DOWNSCAN")
        (currentTime < scanTime + holdTime + flushTime) && (scanState = "UPHOLD")
        (currentTime < scanTime + holdTime) && (scanState = "UPSCAN")
        (currentTime <= holdTime) && (scanState = "HOLD")
        scanState = (dmaState.value == :SMPS) ? scanState : "CLASSIFIER"

        return scanState
    end

    # Set SMPS voltage
    function smps_voltage(t)
        holdTime, scanTime, flushTime, scanLength, startVoltage, endVoltage, c =
            scan_parameters()

        (smps_scan_state.value == "HOLD") && (myV = startVoltage)
        (smps_scan_state.value == "UPSCAN") &&
            (myV = exp(trianglewave(t-holdTime, 2*scanTime, log(startVoltage), log(endVoltage))))
        (smps_scan_state.value == "UPHOLD") && (myV = endVoltage)
        (smps_scan_state.value == "DOWNSCAN") &&
            (myV = exp(trianglewave(t-holdTime-flushTime, 2*scanTime, log(startVoltage), log(endVoltage))))
        (smps_scan_state.value == "FLUSH") && (myV = startVoltage)
        (smps_scan_state.value == "DONE") && (myV = startVoltage)
        (smps_scan_state.value == "CLASSIFIER") && (myV = classifierV.value)

        return myV
    end

    # Generate signals and connect with functions
    smps_scan_state = map(state, smps_elapsed_time)
    V = map(smps_voltage, smps_elapsed_time)
    Dp = map(v -> vtod(Λ, v), V)
    reset = map(
        _ -> push!(smps_start_time, datetime2unix(now(UTC))),
        filter(t -> t > scanLength, smps_elapsed_time),
    )
    smps_scan_state, reset, V, Dp
end

function scan_parameters()
    holdTime = thold.value
    scanTime = tscan.value
    flushTime = tflush.value
    startVoltage = Vhi.value
    endVoltage = Vlow.value
    scanLength = holdTime + 2*scanTime + flushTime
    c = log(endVoltage / startVoltage) / (scanTime)

    return holdTime, scanTime, flushTime, scanLength, startVoltage, endVoltage, c
end
