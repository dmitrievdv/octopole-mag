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

# ax = Axis(fig[1,1])

n_Rin = length(R_ins)
n_W = length(Ws)
n_ang = length(angs)

pos_axs = [Axis(fig[i_Rin, i_ang]) for i_ang in 1:n_ang, i_Rin in 1:n_Rin]
neg_axs = [Axis(fig[i_Rin + n_Rin, i_ang]) for i_ang in 1:n_ang, i_Rin in 1:n_Rin]

model_name = modelname(Ṁ, T_max, R_in, W)

for i_Rin in 1:n_Rin, i_ang in 1:n_ang
    R_in = R_ins[i_Rin]
    ang = angs[i_ang]
    model_name = modelname(Ṁ, T_max, R_in, W)
    pos_ax = pos_axs[i_ang, i_Rin]
    neg_ax = neg_axs[i_ang, i_Rin]
    for Γ in Γs
        Γ_string = rstring(Γ)
        mag = loadmodel(star, "posoct-$(Γ_string)_$model_name")
        prof = HydrogenProfile(star, mag, "Ha_$ang")
        lines!(pos_ax, get_velocities(prof), prof.profile, label = "$Γ")

        mag = loadmodel(star, "negoct-$(Γ_string)_$model_name")
        prof = HydrogenProfile(star, mag, "Ha_$ang")
        lines!(neg_ax, get_velocities(prof), prof.profile, label = "$(-Γ)")
    end
    mag = loadmodel(star, "dipole_$model_name")
    prof = HydrogenProfile(star, mag, "Ha_$ang")
    lines!(pos_ax, get_velocities(prof), prof.profile, color = :red, label = "0.0")
    lines!(neg_ax, get_velocities(prof), prof.profile, color = :red, label = "0.0")
end

axislegend(pos_axs[1], "Positive Γ")
axislegend(neg_axs[1], "Negative Γ")

fig
end