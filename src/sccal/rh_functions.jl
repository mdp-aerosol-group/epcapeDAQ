using DataFrames
using CSV

function ammsulfgf(Dd)
    #   Dd: Dry diameter
    #   Loads Data from E-AIM and converts aw/xw pair to RH/gf pair
    #   RH100 and gf100 are returned as reference for QA purposes
    #   Dd: dry diameter (m)
    #   aw: water activity (-)
    #   k: kappa (-)
    #   vfw: volume fraction of water (-)
    #   gf: hygroscopic growth factor (-)
    #   RH: relative humidity (%)
    #   Note: the Kelvin effect assumes pure water at 298K. Replace
    #   2.1e-9 (m) with "A" parameter from Kohler theory to get T-dependence
    
    data = CSV.read("AIM_aw_xw_gf100nm_trueRH100nm.dat", DataFrame; header = false)
    aw = data[:,1]    
    xw = data[:,2]
    gf100 = data[:,3]
    rh100 = data[:,4]
    asrho = 1769.0 
    asMW = 0.132
    asnu = 3
    h2orho = 997.1
    h2oMW = 0.018
    h2onu = 1
    
    f1 = asrho .* asnu .* h2oMW .* xw
    f2 = h2orho .* asMW .* (1.0 .- xw)
    vfw = f1 ./ (f1 .+ f2)
    kappa = ((aw .- 1.0) .* vfw) ./ ((vfw .- 1) .* aw)
    gf = (1 .+ kappa .* aw ./ (1 .- aw)) .^ (1.0/3.0)
    RH = aw .* exp.(2.1e-9 ./ (gf .* Dd)) * 100.0
    DataFrame(
        aw = aw,
        xw = xw,
        k = kappa,
        gf100 = gf100,
        rh100 = rh100,
        gf = gf,
        RH = RH
    )
end

function getrh(mygf, Dd)
    df = ammsulfgf(Dd)
    ii = argmin((df[!,:gf].-mygf).^2)
    
    gf = df[ii,:gf]
    RH = df[ii,:RH]
    return gf, RH
end

function get_sc(Dd)
    df = ammsulfgf(Dd)
    return maximum(df[!,:RH]) - 100
end

Sc(Dd, κ) = ((4.0 * (2.1e-9)^3.0) / (27.0 .* Dd .^ 3 .* κ))^0.5
sc(Dd, κ) = (exp(Sc(Dd, κ)) - 1) * 100