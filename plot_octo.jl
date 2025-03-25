using TTauUtils
using Makie
using GLMakie

include("names.jl")
include("parameters.jl")

function get_velocities(profile :: HydrogenProfile)
    line_freq = TTauUtils.linefrequency(profile.upper_level, profile.lower_level)

    return -3e5*(profile.frequencies/line_freq .- 1)
end

begin 
fig = Figure()

ax = Axis(fig[1,1])

model_name = modelname(Ṁ, T_max, R_in, W)

for Γ in Γs
    Γ_string = rstring(Γ)
    mag = loadmodel(star, "posoct-$(Γ_string)_$model_name")
    prof = HydrogenProfile(star, mag, "Ha_15")
    lines!(ax, get_velocities(prof), prof.profile)
end

fig
end