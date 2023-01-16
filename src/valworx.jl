U3HANDLE.configIO(EnableCounter0 = false, EnableCounter1 = false, NumberOfTimersEnabled = 0, FIOAnalog = 0)

function valve(pos)
    if pos == true
        U3HANDLE.getFeedback(u3.BitStateWrite(8, 0))
        U3HANDLE.getFeedback(u3.BitStateWrite(9, 1))
    else
        U3HANDLE.getFeedback(u3.BitStateWrite(8, 1))
        U3HANDLE.getFeedback(u3.BitStateWrite(9, 0))
    end
end