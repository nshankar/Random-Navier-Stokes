# Select Fourier coefficients for the vorticity qhat in the domain [0, N//2] x [-N//2, N//2] and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles
using FFTW
using Plots

function gauss(i,j,center, var, coef)
	return coef * exp(-((center[1] - i)^2 + (center[2] - j)^2)/(2*var))

end

filename = "initialConditions/vorticityFreqIC5.csv"

centers = [[6,6], [20, 30], [44, 15]]
coefs = [3, 5, -4]
vars = [1, 10, 4]

N = 61 # Require N odd (for now)
q(i,j) = gauss(i,j, centers[1], vars[1], coefs[1]) + gauss(i,j, centers[2], vars[2],coefs[2]) + gauss(i,j, centers[3], vars[3], coefs[3])

ICReal = Array{Float64}(undef,(N,N))
for i=1:N
	for j=1:N
		ICReal[i,j] = q(i,j)
	end
end

IC = rfft(ICReal)

#=
maxFreq = Int((N-1)/2)
IC = Array{ComplexF64}(undef, (Int((N+1)/2), N))

for i=0:maxFreq
	for j=-maxFreq:maxFreq
		IC[i + 1, mod(j, N) + 1] = qhat(i, j)
	end
end
IC[1,1] = 0
=#
heatmap(real.(IC), show=true)

writedlm(filename,  IC, ',')

