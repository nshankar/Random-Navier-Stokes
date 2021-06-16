# Select Fourier coefficients for the vorticity qhat in the domain [0, N//2] x [-N//2, N//2] and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles

A = Set([1,5,6,8,17])

N = 31 # Require N odd (for now)
qhat(i,j) = 100*( (i in A && j in A) || mod(i*j, N)+1 in A)
filename = "vorticityFreqIC2.csv"

maxFreq = Int((N-1)/2)
IC = Array{Complex{Float64}}(undef, (Int((N+1)/2), N))

for i=0:maxFreq
	for j=-maxFreq:maxFreq
		IC[i + 1, mod(j, N) + 1] = qhat(i, j)
	end
end

writedlm(filename,  IC, ',')