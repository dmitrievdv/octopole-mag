using TTauUtils
using Makie
using GLMakie
using LeastSquaresOptim

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

function get_inc_fitting_function(prof_to_fit)
    function fitting_function(pars)
        println(pars)
        mag = StationarySolidMagnetosphereNHCool("fit_mag", star, pars[2], pars[2] + pars[3], 10^(pars[4]), T_max, 10; mag_kwargs...)
        prof = HydrogenProfile(mag, 3, 2, pars[1], prof_precision_args...; prof_kwargs...)
        return prof.profile .- prof_to_fit.profile
    end
    return fitting_function
end

res = optimize(get_inc_fitting_function(prof_to_fit), [inc_ang*1.1, 4.0*1.1, 2.0*1.1, -8.0], LevenbergMarquardt())
pars = res.minimizer

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