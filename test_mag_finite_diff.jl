using TTauUtils
using Makie
using GLMakie

include("names.jl")
include("parameters.jl")

R_in = 4.0
W = 2.0
Ṁ = 1e-8
T_max = 7000 

R_in = 4
W = 2
Ṁ = 1e-8
T_max = 8000
inc_ang = 75

prof_precision_args = (0.02, 0.05, 0.1, 30)
prof_kwargs = (:voigt => false, :progress_output => false, :blue_v_max => 300)
mag_kwargs = (:progress_output => false, n_t = 10, n_r_m = 4)

mag = StationarySolidMagnetosphereNHCool("mag", star, R_in, R_in + W, Ṁ, T_max, 10; mag_kwargs...)
prof = HydrogenProfile(mag, 3, 2, inc_ang, prof_precision_args...; prof_kwargs...)

naught_val = 0.0

eps = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6]
n_eps = length(eps)
finite_diffs = zeros(n_eps)

for i_eps = 1:n_eps
    println(eps[i_eps])
    shifted_mag = StationarySolidMagnetosphereNHCool("mag_shifted", star, R_in, R_in + W, Ṁ, T_max*(1+eps[i_eps]), 10; mag_kwargs...)
    shifted_prof = HydrogenProfile(shifted_mag, 3, 2, inc_ang, prof_precision_args...; prof_kwargs...)
    shifted_val = sum((shifted_prof.profile .- prof.profile) .^ 2)
    finite_diffs[i_eps] = shifted_val/eps[i_eps]
end

fig = Figure()
ax = Axis(fig[1,1], xscale = log10, yscale = log10, title = "T_max")
scatter!(ax, eps, finite_diffs)

fig
