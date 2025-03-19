function rstring(r; max_d = 2)
    dec_char = ["d", "c"]
    r_int = max_d
    while (round(Int, r*10^(r_int*1.0)) % 10) == 0
        r_int -= 1
    end
    r_int
    whole = floor(Int, r)
    dec_str = if r_int <= 0
        string(whole)
    else 
        part = floor(Int, r*10^r_int) - whole*10^r_int
        zeros = r_int - floor(Int, log10(part)) - 1
        dec_char[r_int]*string(whole)*("0"^zeros)*string(part)
    end
    return dec_str
end

function modelname(Ṁ, T_max, R_in, W)
    mag_string = "$(rstring(R_in))-$(rstring(R_in+W))"
    T_string = string(round(Int, T_max))
    lg_Ṁ = log10(Ṁ)
    Ṁ_string = string(round(Int, -lg_Ṁ*10))
    return "$(Ṁ_string)_$(T_string)_$(mag_string)"
end

function linename(u, l)
    line_name = ""
    if l == 2
        line_name *= 'H'
    end
    if l == 3
        line_name *= "Pa"
    end
    if l == 4
        line_name *= "Br"
    end
    line_name *= 'a' + u-1-l
end