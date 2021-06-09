# Select Fourier coefficients for the vorticity qhat in the domain [0, N//2] x [-N//2, N//2] and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles

qhat(i,j) = (i == 1) && (j == 0 || j == 1)
N = 11 # Require N odd (for now)
filename = "vorticityFreqIC.csv"

maxFreq = Int((N-1)/2)
IC = Array{Float64}(undef, (Int((N+1)/2), N))

for i=0:maxFreq
	for j=-maxFreq:maxFreq
		IC[i + 1, mod(j, N) + 1] = qhat(i, j)
	end
end

writedlm(filename,  IC, ',')