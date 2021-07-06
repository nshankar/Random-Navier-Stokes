# Select Fourier coefficients for the vorticity qhat in the domain [0, N//2] x [-N//2, N//2] and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles

A = Set([2, 3, 5, 8, 13, 21])

N = 31 # Require N odd (for now)
qhat(i,j) = 10*((i == 2 || i == 3) && (j == 5))

#10*( (i in A && j in A) || mod(i*j, N)+1 in A)

filename = "initialConditions/vorticityFreqIC3.csv"

maxFreq = Int((N-1)/2)
IC = Array{ComplexF64}(undef, (Int((N+1)/2), N))

for i=0:maxFreq
	for j=-maxFreq:maxFreq
		IC[i + 1, mod(j, N) + 1] = qhat(i, j)
	end
end
IC[1,1] = 0

writedlm(filename,  IC, ',')