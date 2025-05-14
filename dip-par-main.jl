using TTauUtils
using Makie
using GLMakie

include("names.jl")
include("parameters.jl")

function get_velocities(profile :: HydrogenProfile)
    line_freq = TTauUtils.linefrequency(profile.upper_level, profile.lower_level)

    return -3e5*(profile.frequencies/line_freq .- 1)
end


Γ = 3

prof_precision_args = (0.02, 0.05, 0.1, 100)
prof_kwargs = (:voigt => true, :progress_output => false)

mag_kwargs = (:progress_output => false,)

model_name = modelname(Ṁ, T_max, R_in, W)
Γ_string = rstring(Γ)
octo_mag = StationaryOctoDipoleMagnetosphereNHCool("posoct-$(Γ_string)_$model_name", star, Ṁ, T_max, R_in, R_in + W, Γ, 10; mag_kwargs...)
inc_ang = 45
octo_prof = HydrogenProfile(octo_mag, 3, 2, inc_ang, prof_precision_args...; prof_kwargs...)

n_T_max_dip = length(dip_T_maxs)
n_W_dip = length(dip_Ws)
n_R_in_dip = length(dip_R_ins)
n_inc_ang_dip = length(dip_angs)

# parameter_shape = (n_T_max_dip, n_R_in_dip, n_W_dip, n_inc_ang_dip)
parameter_arrays = (dip_T_maxs, dip_R_ins, dip_Ws, dip_angs)
parameter_shape = Tuple(length(arr) for arr in parameter_arrays)
n_parameters = length(parameter_shape)

parameter_array = zeros(parameter_shape..., n_parameters)


for parameter_id in 1:n_parameters
    for i_parameter = 1:parameter_shape[parameter_id]
        selectdim(selectdim(parameter_array, parameter_id, i_parameter), n_parameters, parameter_id) .= parameter_arrays[parameter_id][i_parameter]
    end
end


δ_array = zeros(parameter_shape...)
parameter_array[5, 2, 3, 15, :]

n_freqs = prof_precision_args[end]
dip_profile_array = zeros(parameter_shape..., n_freqs)
dip_velocities_array = zeros(parameter_shape..., n_freqs)

for parameter_index in CartesianIndices(δ_array)
    T_max, R_in, W, inc_ang = parameter_array[parameter_index, :]
    dip_model_name = "dip_"*modelname(Ṁ, T_max, R_in, W)
    dip_mag = loadmodel(star, dip_model_name)
    dip_prof_name = "Ha$(round(Int, inc_ang))"
    dip_prof = HydrogenProfile(star, dip_mag, dip_prof_name)
    dip_profile_array[parameter_index, :] .= dip_prof.profile
    dip_velocities_array[parameter_index, :] .= get_velocities(dip_prof)
    δ_array[parameter_index] = sqrt(sum((dip_prof.profile - octo_prof.profile).^2)/n_freqs)
end

min_δ, min_index = findmin(identity, δ_array)
dip_parameter_index = CartesianIndex(6, 6, 6, 12)

fig = Figure()
ax = Axis(fig[1,1])

lines!(ax, get_velocities(octo_prof), octo_prof.profile, label = "octo")
lines!(ax, dip_velocities_array[dip_parameter_index, :], dip_profile_array[dip_parameter_index, :], label = "dip")
lines!(ax, dip_velocities_array[min_index, :], dip_profile_array[min_index, :], label = "dip min")
axislegend(ax)

fig

good_pars = [parameter_array[index, :] for index in findall(x -> x < √2*min_δ, δ_array)]
mean_pars = sum(good_pars)/length(good_pars)

# calculate_δ!(δ_array, parameter_array, parameter_arrays, )
