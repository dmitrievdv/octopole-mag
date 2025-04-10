using TTauUtils

include("names.jl")
include("parameters.jl")

prof_precision_args = (0.02, 0.05, 0.1, 100)
prof_kwargs = (:voigt => true, :progress_output => false)

mag_kwargs = (:progress_output => false,)

for R_in in R_ins, W in Ws
    model_name = modelname(Ṁ, T_max, R_in, W)
    dip_mag = StationarySolidMagnetosphereNHCool("dipole_$model_name", star, R_in, R_in + W, Ṁ, T_max, 10; mag_kwargs...)
    savemodel(dip_mag)  
    for i_ang in angs
        prof_name = linename(u, l)*"_$i_ang"
        println(prof_name)
        prof = HydrogenProfile(dip_mag, u, l, i_ang, prof_precision_args...; prof_kwargs...)
        saveprofile(prof, prof_name)
    end

    for Γ in Γs
        println(Γ)
        Γ_string = rstring(Γ)
        octo_mag = StationaryOctoDipoleMagnetosphereNHCool("posoct-$(Γ_string)_$model_name", star, Ṁ, T_max, R_in, R_in + W, Γ, 10; mag_kwargs...)
        savemodel(octo_mag)
    
        for i_ang in angs
            prof_name = linename(u, l)*"_$i_ang"
            println(prof_name)
            prof = HydrogenProfile(octo_mag, u, l, i_ang, prof_precision_args...; prof_kwargs...)
            saveprofile(prof, prof_name)
        end

        println(-Γ)
        octo_mag = StationaryOctoDipoleMagnetosphereNHCool("negoct-$(Γ_string)_$model_name", star, Ṁ, T_max, R_in, R_in + W, -Γ, 10; mag_kwargs...)
        savemodel(octo_mag)
        for i_ang in angs
            prof_name = linename(u, l)*"_$i_ang"
            println(prof_name)
            prof = HydrogenProfile(octo_mag, u, l, i_ang, prof_precision_args...; prof_kwargs...)
            saveprofile(prof, prof_name)
        end
    end
end
