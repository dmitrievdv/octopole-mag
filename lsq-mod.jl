using TTauUtils
using Makie
using GLMakie
import LeastSquaresOptim as lsqopt
import Optim as opt

include("names.jl")
include("parameters.jl")

function get_velocities(profile :: HydrogenProfile)
    line_freq = TTauUtils.linefrequency(profile.upper_level, profile.lower_level)

    return -3e5*(profile.frequencies/line_freq .- 1)
end

Γ = 0.01
R_in = 4
W = 2
Ṁ = 1e-8
T_max = 7000
inc_ang = 75

prof_precision_args = (0.02, 0.05, 0.1, 30)
prof_kwargs = (:voigt => false, :progress_output => false, :blue_v_max => 300)
mag_kwargs = (:progress_output => false,)

mag_to_fit = StationaryOctoDipoleMagnetosphereNHCool("mag_to_fit", star, Ṁ, T_max, R_in, R_in + W, Γ, 10; mag_kwargs...)
prof_to_fit = HydrogenProfile(mag_to_fit, 3, 2, inc_ang, prof_precision_args...; prof_kwargs...)

function get_inc_fitting_function(prof_to_fit, pars_start, pars_scale)
    function fitting_function(pars)
        pars_real = pars .* pars_scale .+ pars_start
        println(pars, " ", pars_real)
        mag = StationarySolidMagnetosphereNHCool("fit_mag", star, pars_real[2], pars_real[2] + pars_real[3], 10^(pars_real[4]), T_max, 10; mag_kwargs...)
        prof = HydrogenProfile(mag, 3, 2, pars[1], prof_precision_args...; prof_kwargs...)
        return prof.profile .- prof_to_fit.profile
    end
    return fitting_function
end

function get_fitting_function_grad(prof_to_fit, pars_start, pars_scale, eps = 0.01)
    function fitting_function_grad!(grad, pars)
        n_pars = length(pars)
        pars_real = pars .* pars_scale .+ pars_start
        println(pars, " ", pars_real)
        mag = StationarySolidMagnetosphereNHCool("fit_mag", star, pars_real[2], pars_real[2] + pars_real[3], 10^(pars_real[4]), T_max, 10; mag_kwargs...)
        prof = HydrogenProfile(mag, 3, 2, pars[1], prof_precision_args...; prof_kwargs...)

        naught_value = sum((prof.profile .- prof_to_fit.profile) .^ 2)

        pars_real_shifted = (pars .+ eps) .* pars_scale .+ pars_start
        pars_dim = zeros(n_pars)
        pars_dim .= pars_real
        for i_par = 1:n_pars
            pars_dim[i_par] = pars_real_shifted[i_dim]
            mag = StationarySolidMagnetosphereNHCool("fit_mag", star, pars_dim[2], pars_dim[2] + pars_dim[3], 10^(pars_dim[4]), T_max, 10; mag_kwargs...)
            prof = HydrogenProfile(mag, 3, 2, pars_dim[1], prof_precision_args...; prof_kwargs...)
            pars_dim[i_par] = pars_real[i_par]
            # shifted_values[i_par] = sum((prof.profile .- prof_to_fit.profile) .^ 2)
            grad[i_par] = (sum((prof.profile .- prof_to_fit.profile) .^ 2) - naught_value)/eps
        end
    end
    return fitting_function_grad
end

pars_start = [inc_ang*0.9, 4.0*0.9, 2.0*0.9, -8.0]
pars_scale = [10.0, 1.0, 1.0, 0.1]
func = get_inc_fitting_function(prof_to_fit, pars_start, pars_scale)
grad = get_fitting_function_grad(prof_to_fit, pars_start, pars_scale, 0.01)
res = opt.optimize(v -> sum( func(v) .^ 2), grad, zeros(length(pars_start)))
pars_res = opt.minimizer(res)

pars = pars_res .* pars_scale .+ pars_start
fitted_inc = pars[1]

dip_mag = StationarySolidMagnetosphereNHCool("dip_mag", star, pars[2], pars[2] + pars[3], 10^(pars[4]), T_max, 10; mag_kwargs...)
dip_prof = HydrogenProfile(dip_mag, 3, 2, fitted_inc, prof_precision_args...; prof_kwargs...)
# dip_prof_same = HydrogenProfile(dip_mag, 3, 2, fitted_inc*1.1, prof_precision_args...; prof_kwargs...)

fig = Figure()
ax = Axis(fig[1,1])

lines!(ax, get_velocities(prof_to_fit), prof_to_fit.profile)
lines!(ax, get_velocities(dip_prof), dip_prof.profile)
lines!(ax, get_velocities(dip_prof_same), dip_prof_same.profile)
fig