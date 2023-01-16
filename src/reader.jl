using DelimitedFiles
using Chain
using Plots
using Dates

s = open("all.txt") do file
    read(file, String)
end

b = split(s, "\r")
txcpc = split(s, "\n")
tcpc = filter(x -> length(x) > 1, txcpc)
tccn = filter(x -> String(x[1:2]) == ";H", b)
tpops = filter(x -> String(x[1:2]) == ";P", b)
txlj = filter(x -> String(x[1:2]) == ";L", b)
tlj = map(x -> split(x, "\n")[1], txlj)

t = map(x -> (@chain split(x, ";") getindex(_, 1) DateTime), tcpc) 
cc = map(x -> (@chain split(x, ",") getindex(_, 9) parse(Float64, _)), tlj) 
Vr = map(x -> (@chain split(x, ",") getindex(_, 5) parse(Float64, _)), tlj) 
ccn = map(x -> (@chain split(x, ",") getindex(_, 21) parse(Float64, _)), tccn) 

tm = t[1]:Minute(2):t[end]
ticks = Dates.format.(tm, "HH:MM")

ii = 1:1000
p1 = plot(t[ii], cc[ii], xticks=(tm,ticks), label = "CN")
p1 = plot!(t[ii], ccn[ii], label = "CCN")
p2 = plot(t[ii], Vr[ii],xticks=(tm,ticks), label = "V")

plot(p1, p2, layout = grid(2,1))