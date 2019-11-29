using LinearAlgebra
using SkyCoords

ra = 89.014303
dec = 13.924912
ω = 37.59

α = deg2rad(ra)
δ = deg2rad(dec)
d = 1000.0 / ω
println("α=",α,"\tδ=",δ,"\td=",d)

pmra = 372.72
pmdec = -483.69
radial_velocity = 0.37

vSUN = [11.1, 232.24, 7.25]#km/s

tropical_year = 31556925.445#seconds in one tropical year
parsec = 30856775812800#km in one parsec
factor = parsec / tropical_year
println("parsec,year conversion factor = ", factor)

#SI_unit = 2.44507E-17
#println("SI mas/year conversion unit:", SI_unit)

#=
μRA = pmra / cos(δ)
μDEC = pmdec

vRA = 4.74 * d * μRA / 1000.0#km/s
vDEC = 4.74 * d * μDEC / 1000.0#km/s
vR = radial_velocity#km/s

println("vRA:", vRA, "\tvREC:", vDEC, "vR:", vR)
println("adj. vRA:", vRA+vSUN[1], "\tadj. vREC:", vDEC+vSUN[2], "adj. vR:", vR-vSUN[3])
=#

η = deg2rad(58.5986320306)
αGC = deg2rad(266.4051)
δGC = deg2rad(-28.936175)
dGC = 8300
xGC = [1, 0, 0]
z = 27
θ = asin(z / dGC)
println("θ=",θ)

#positions
r_icrs = d * [cos(α) * cos(δ), sin(α) * cos(δ), sin(δ)]

#velocities
#pmra = 4.74 * d * pmra / 1000.0#in line with km/s
#pmdec = 4.74 * d * pmdec / 1000.0#in line with km/s
#v_icrs = d * [cos(α) * cos(δ), sin(α) * cos(δ), sin(δ)]
pmra = pmra / tropical_year / 1000
pmdec = pmdec / tropical_year / 1000
d *= parsec

V_icrs = [radial_velocity*cos(α)*cos(δ)-d*sin(α)*pmra-d*cos(α)*sin(δ)*pmdec, radial_velocity*sin(α)*cos(δ)+d*cos(α)*pmra-d*sin(α)*sin(δ)*pmdec, radial_velocity*sin(δ)+d*cos(δ)*pmdec]
println("V_icrs=",V_icrs)

R₁ = Matrix(undef, 3, 3)
R₁[1,1] = cos(δGC)
R₁[1,2] = 0
R₁[1,3] = sin(δGC)
R₁[2,:] = [0 1 0]
R₁[3,1] = - sin(δGC)
R₁[3,2] = 0
R₁[3,3] = cos(δGC)
println("R₁:", R₁)

R₂ = Matrix(undef, 3, 3)
R₂[1,1] = cos(αGC)
R₂[1,2] = sin(αGC)
R₂[1,3] = 0
R₂[2,1] = - sin(αGC)
R₂[2,2] = cos(αGC)
R₂[2,3] = 0
R₂[3,:] = [0 0 1]
println("R₂:", R₂)

R₃ = Matrix(undef, 3, 3)
R₃[1,:] = [1 0 0]
R₃[2,1] = 0
R₃[2,2] = cos(η)
R₃[2,3] = sin(η)
R₃[3,1] = 0
R₃[3,2] = - sin(η)
R₃[3,3] = cos(η)
println("R₃:", R₃)

R = R₃ * R₁ * R₂
println("R:", R)

R = Matrix(undef, 3, 3)
R[1,1] = cos(αGC) * cos(δGC)
R[1,2] = cos(δGC) * sin(αGC)
R[1,3] = sin(δGC)
R[2,1] = - cos(αGC) * sin(δGC) * sin(η) - sin(αGC) * cos(η)
R[2,2] = - sin(αGC) * sin(δGC) * sin(η) + cos(αGC) * cos(η)
R[2,3] = cos(δGC) * sin(η)
R[3,1] = - cos(αGC) * sin(δGC) * cos(η) + sin(αGC) * sin(η)
R[3,2] = - sin(αGC) * sin(δGC) * cos(η) - cos(αGC) * sin(η)
R[3,3] = cos(δGC) * cos(η)
println("R:", R)

H = Matrix(undef, 3, 3)
H[1,1] = cos(θ)
H[1,2] = 0
H[1,3] = sin(θ)
H[2,1] = 0
H[2,2] = 1
H[2,3] = 0
H[3,1] = - sin(θ)
H[3,2] = 0
H[3,3] = cos(θ)
println("H:", H)
H = [cos(θ) 0 sin(θ); 0 1 0; - sin(θ) 0  cos(θ)]

r = R * r_icrs - dGC * xGC
rGC = H * r

println("distance:", d)
println("r_icrs:", r_icrs)
println("r:", r)
println("rGC:", rGC)

#+++++++++++++ velocities ++++++++++++++++++++++++
c1 = ICRSCoords(α, δ)
c2 = convert(GalCoords, c1)
l = c2.l
b = c2.b
println("c1:", c1, "\tc2:", c2)

#transform angular velocities to the galactic frame
C₁ = sin(δGC) * cos(δ) - cos(δGC) * sin(δ) * cos(α - αGC)
C₂ = cos(δGC) * sin(α - αGC)
cosb = sqrt(C₁^2 + C₂^2)

C = Matrix(undef, 2, 2)
C[1,1] = C₁
C[1,2] = C₂
C[2,1] = - C₂
C[2,2] = C₁

μ = [pmra, pmdec]
μ = C * μ / cosb
println(μ)
μl = μ[1]
μb = μ[2]

Vl = 4.74 * μl * d / 1000.0
Vb = 4.74 * μb * d / 1000.0
Vᵣ = radial_velocity

U = Vᵣ * cos(l) * cos(b) - Vl * sin(l) - Vb * cos(l) * sin(b)
V = Vᵣ * sin(l) * cos(b)+ Vl * cos(l) - Vb * sin(l) * sin(b)
W = Vᵣ * sin(b) + Vb * cos(b) 
println("U:", U, "\tV:", V, "\tW:", W)
#println("adj. U:", U+vSUN[1], "\tadj. V:", V+vSUN[2], "\tadj. W:", W+vSUN[3])

VGC = [28.46188771, 157.93916086, 17.65189313]
#X = H * R * V_icrs - VGC
#X = R * V_icrs - inv(H) * VGC
println("VGC1 = ",VGC)

#X=[3701.88, -15812.7, 2195.97]
#VGC = H * R * V_icrs - X
VGC = H * R * V_icrs
println("VGC2 = ", VGC)