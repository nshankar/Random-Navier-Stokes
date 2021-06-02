# Select a periodic velocity field (u1(x,y), u2(x,y)) on [0,2*pi]^2 and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles

u1(x,y) = exp(-x^2 - y^2) * x * y
u2(x,y) = sin(exp(x*y))
N = 201 # Require N odd (for now)
filename = "IC.csv"


IC = Array{Float64}(undef, (N,N,2))

for i=1:N
	for j=1:N
		IC[i,j,1] = u1(2*pi*i/N, 2*pi*j/N)
		IC[i,j,2] = u2(2*pi*i/N, 2*pi*j/N)
	end
end

writedlm(filename,  IC, ',')