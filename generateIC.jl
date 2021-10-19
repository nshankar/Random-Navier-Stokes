# Select Fourier coefficients for the vorticity qhat in the domain [0, N//2] x [-N//2, N//2] and a resolution N
# Generates a csv file which can be used as input for main.jl
using DelimitedFiles
using FFTW

function gauss(i,j,center, var, coef)
	return coef * exp(-((center[1] - i)^2 + (center[2] - j)^2)/(2*var))

end

function generateIC()
	filename = "initialConditions/vortexPair.csv"


	N = 101 # Require N odd
	#q(i,j) = gauss(i,j, [10, 10], 2^2, -10) + gauss(i,j, [65, 40], 10^2, -5) + gauss(i,j, [20, 80], 5^2, 12)

	a = 35.5
	b = 14.5
	var = 6^2
	q(i,j) = gauss(i,j, [a, b], var, -5) + gauss(i,j, [N-a, b], var, 5) + 
		+ gauss(i,j, [a, N-b], var, 5) + gauss(i,j, [N-a, N-b], var, -5)

	ICReal = Array{Float64}(undef,(N,N))

	for i=1:N
		for j=1:N
			ICReal[i,j] = q(i,j)
		end
	end


	IC = rfft(ICReal)
	IC[1,1] = 0

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

	writedlm(filename,  IC, ',')
end

generateIC()

