using TTauUtils
using Makie
using GLMakie

include("parameters.jl")
include("names.jl")
# include("kwargs.jl")

prof_precision_args = (0.02, 0.05, 0.1, 100)
prof_kwargs = (:voigt => false, :progress_output => true, blue_v_max = 300)

mag_kwargs = (:progress_output => false,)

calc_models = true

lines = [(5,2), (6,2), (7,2)] 
n_lines = length(lines)

T_max = 8000

model_name = modelname(Ṁ, T_max, R_in, W)


dip_mag = if calc_models 
    StationarySolidMagnetosphereNHCool("dipole_"*model_name, star, R_in, R_in + W, Ṁ, T_max, 10; mag_kwargs...)
else
    dip_mag = loadmodel(star, "dipole_"*model_name)
end

if calc_models
    savemodel(dip_mag)
end



ang = 15

for line in lines
    println("Dipole profile...")
    dip_prof = HydrogenProfile(dip_mag, line..., ang, prof_precision_args...; prof_kwargs...)
    line_string = linename(line...)
    saveprofile(dip_prof, "$(line_string)_$ang")
end

Γ = 2

posoct_mag = if calc_models
    StationaryOctoDipoleMagnetosphereNHCool("posoct-$(rstring(abs(Γ)))_"*model_name, star, Ṁ, T_max, R_in, R_in + W, Γ, 10; mag_kwargs...)
else
    loadmodel(star, "posoct-$(rstring(abs(Γ)))_"*model_name)
end

negoct_mag = if calc_models 
    StationaryOctoDipoleMagnetosphereNHCool("negoct-$(rstring(abs(Γ)))_"*model_name, star, Ṁ, T_max, R_in, R_in + W, -Γ, 10; mag_kwargs...)
else
    loadmodel(star, "negoct-$(rstring(abs(Γ)))_"*model_name)
end

if calc_models
    savemodel(posoct_mag)
    savemodel(negoct_mag)
end


for line in lines
    println("Posoct profile...")
    posoct_prof = HydrogenProfile(posoct_mag, line..., ang, prof_precision_args...; prof_kwargs...)
    line_string = linename(line...)
    saveprofile(posoct_prof, "$(line_string)_$ang")

    println("Negoct profile...")
    negoct_prof = HydrogenProfile(negoct_mag, line..., ang, prof_precision_args...; prof_kwargs...)
    line_string = linename(line...)
    saveprofile(negoct_prof, "$(line_string)_$ang")
end

function get_velocities(profile :: HydrogenProfile)
    line_freq = TTauUtils.linefrequency(profile.upper_level, profile.lower_level)
    return -3e5*(profile.frequencies/line_freq .- 1)
end

fig = Figure()


axs = [Axis(fig[1,i_line]) for i_line = 1:n_lines]

prof = HydrogenProfile(star, dip_mag, "$(linename(lines[1]...))_$ang")
vels = get_velocities(prof)

for (i_line, line) in enumerate(lines)
    dip_prof = HydrogenProfile(star, dip_mag, "$(linename(line...))_$ang")
    posoct_prof = HydrogenProfile(star, posoct_mag, "$(linename(line...))_$ang")
    negoct_prof = HydrogenProfile(star, negoct_mag, "$(linename(line...))_$ang")

    lines!(axs[i_line], vels, dip_prof.profile, color = :black)
    lines!(axs[i_line], vels, posoct_prof.profile, color = :red)
    lines!(axs[i_line], vels, negoct_prof.profile, color = :blue)
end

fig