using Distributed
import REPL

terminal = REPL.Terminals.TTYTerminal(string(), stdin, stdout, stderr)

n_proc = 4
addprocs(n_proc)

@everywhere using TTauUtils


@everywhere include("names.jl")
@everywhere include("parameters.jl")

prof_precision_args = (0.02, 0.05, 0.1, 100)
prof_kwargs = (:voigt => true, :progress_output => false)

mag_kwargs = (:progress_output => false,)

clear_console() = print("\033c")

n_mags = length(dip_R_ins)*length(dip_Ws)
i_mag = 1

for R_in in dip_R_ins, W in dip_Ws
    println("Calculating models... $i_mag of $n_mags : R_in = $R_in, W = $W")
    mags = @distributed vcat for T_max in dip_T_maxs
        mag = StationarySolidMagnetosphereNHCool("dip_"*modelname(Ṁ, T_max, R_in, W), star, R_in, R_in + W, Ṁ, T_max, 10; mag_kwargs...)
    end
    println("Saving models...")
    for mag in mags
        savemodel(mag)
    end

    global i_mag += 1
    # mag_percent = round(Int, i_mag*100/n_mags)

    i_profs = 1

    for T_max in dip_T_maxs
        println("Calculating profiles... $i_profs of $(length(dip_T_maxs)) : R_in = $R_in, W = $W, T_max = $T_max")
        mag = loadmodel(star, "dip_"*modelname(Ṁ, T_max, R_in, W))
        profs = @distributed vcat for ang in dip_angs
            prof = HydrogenProfile(mag, u, l, ang, prof_precision_args...; prof_kwargs...)
        end
        i_profs += 1
        println("Saving profiles...")
        for prof in profs
            ang = round(Int, prof.orientation.i / π * 180)
            prof_name = linename(u, l)*"$ang"
            saveprofile(prof, prof_name)
        end
        REPL.Terminals.cmove_line_up(terminal); REPL.Terminals.clear_line(terminal)
        REPL.Terminals.cmove_line_up(terminal); REPL.Terminals.clear_line(terminal)
    end
    REPL.Terminals.cmove_line_up(terminal); REPL.Terminals.clear_line(terminal)
    REPL.Terminals.cmove_line_up(terminal); REPL.Terminals.clear_line(terminal)
end

