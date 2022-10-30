using Plots
using DataFrames
using Interpolations
using CSV

function get_points(V)
    println(V)
    push!(classifierV, V)

    sleep(10)
    x = map(1:10) do i
        sleep(1)
        labjack_signals.value[1][3] |> (x -> x * 1000)
    end

    return DataFrame(setV=V, readV=x)
end

df = mapfoldl(get_points, vcat, [10:10:90; 100:100:1000; 2000:1000:10000])
df1 = filter(:readV => x -> x > 10, df)
df2 = sort(df1, [:readV])

itp = interpolate((df2[!, :readV],), df2[!, :setV], Gridded(Linear()))
xdata = 10:10.0:10000.0
extp = extrapolate(itp, Flat())
p = scatter(df1[!, :readV],
    df1[!, :setV], 
    xscale=:log10, 
    yscale=:log10,
    xlim=(10, 10000), 
    ylim=(10, 10000),
    xlabel = "Read Voltage (V)",
    ylabel = "Set Voltage (V)",
    legend = :bottomright,
    color = :darkred,
    label = "Data"
)

p = plot!(xdata, extp.(xdata), color = :black, label = "Fit")

df2 |>  CSV.write("voltage_calibration.csv")
savefig(p, "voltage_calibration.pdf")
display(p)