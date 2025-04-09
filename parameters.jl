Ṁ = 1e-8
T_max = 7000

u = 3; l = 2

R_in = 4.0
W = 2.0

R_ins = [4]#[4, 5, 6]
Ws = [0.5, 1, 2]
Γs = [0.5, 1, 2, 4, 8]

max_Γ = 2/3*R_in^2

angs = [15, 35, 55, 75]

dip_R_ins = [3:0.2:5;]
dip_Ws = [1:0.2:3;]
dip_T_maxs = [6500:100:7500;]
dip_angs = [1:4:89;]

star = Star("cTTS", 2, 0.5, 4000, 10)
savestar(star)