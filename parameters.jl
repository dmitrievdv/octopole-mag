Ṁ = 1e-8
T_max = 7000

u = 3; l = 2

R_in = 4.0
W = 2.0

model_name = modelname(Ṁ, T_max, R_in, W)

max_Γ = 2/3*R_in^2

angs = [45]

star = Star("cTTS", 2, 0.5, 4000, 10)
savestar(star)