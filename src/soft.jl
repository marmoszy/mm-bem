using LinearAlgebra
#using IterativeSolvers

# Input parameters
fname = length(ARGS)>0 ? ARGS[1] : "sphere-1.905-600.msh"
th = length(ARGS)>1 ? parse(Float32,ARGS[2]) : 0.0
f0 = length(ARGS)>2 ? parse(Float32,ARGS[3]) : 38e3
c0 = length(ARGS)>3 ? parse(Float32,ARGS[4]) : 1480.0
k = 2π * f0 / c0
d = [cos(π/180 * th), sin(π/180 * th), 0.0]

# Read mesh file
lines = readlines(fname)
n = parse(Int, lines[5])
v = [parse.(Float64, split(lines[5 + i])[2:4]) for i in 1:n]
m = parse(Int, lines[n + 8])
e = [parse.(Int, split(lines[n + 8 + i])[6:8]) for i in 1:m]

# Triangle centers (P0)
x = [sum([v[e[i][j]] / 3 for j in 1:3]) for i in 1:m]
# Three points from each triangle (P1)
# m *=3; x = [sum([v[e[i÷3+1][j+1]]*(i%3==j ? 4 : 1)/6 for j in 0:2]) for i in 0:m-1]

# Helmholtz single layer operator
S = zeros(ComplexF64, m, m)
for i in 1:m
    for j in 1:m
        r = norm(x[i] - x[j])
        S[i, j] = r != 0 ? exp(im * k * r) / r / (4π) : im * k / (4π)
    end
end

# Incident field
f = [exp(im * k * x[i]' * d) for i in 1:m]

# Surface solution
φ = lu(S) \ (-f)
#φ = gmres(S,-f)

# Far-field directions
r1 = [[cos(π/180 * i), sin(π/180 * i), 0.0] for i in 0:359]
psc = [sum([exp(-im * k * x[j]' * r) / (4π) * φ[j] for j in 1:m]) for r in r1]

# Print scattering pattern
for i in 1:360
    println(i - 1, '\t', abs(psc[i]))
end
