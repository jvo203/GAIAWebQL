ra = 89.014303
dec = 13.924912
ω = 37.59

α = deg2rad(ra)
δ = deg2rad(dec)
d = 1000.0 / ω

pmra = 372.72
pmdec = -483.69
radial_velocity = 0.37

η = deg2rad(58.5986320306)
αGC = deg2rad(266.4051)
δGC = - deg2rad(-28.936175)#the sign in the rotation matrix as per documentation seems to be inverted; cannot understand the root cause of the problem
dGC = 8300
xGC = [1, 0, 0]
z = 27
θ = asin(z / dGC)

r_icrs = d * [cos(α) * cos(δ), sin(α) * cos(δ), sin(δ)]

R₁ = Matrix(undef, 3, 3)
R₁[1,1] = cos(δGC)
R₁[1,2] = 0
R₁[1,3] = - sin(δGC)
R₁[2,:] = [0 1 0]
R₁[3,1] = sin(δGC)
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
R[1,3] = - sin(δGC)
R[2,1] = cos(αGC) * sin(δGC) * sin(η) - sin(αGC) * cos(η)
R[2,2] = sin(αGC) * sin(δGC) * sin(η) + cos(αGC) * cos(η)
R[2,3] = cos(δGC) * sin(η)
R[3,1] = cos(αGC) * sin(δGC) * cos(η) + sin(αGC) * sin(η)
R[3,2] = sin(αGC) * sin(δGC) * cos(η) - cos(αGC) * sin(η)
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

r = R * r_icrs - dGC * xGC
rGC = H * r

println("distance:", d)
println("r_icrs:", r_icrs)
println("r:", r)
println("rGC:", rGC)
