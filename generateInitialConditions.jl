using DelimitedFiles
using FFTW

function gauss(i,j,center, var, coef)
	return coef * exp(-((center[1] - i)^2 + (center[2] - j)^2)/(2*var))

end

"""
Modify as needed to generate a csv file which can be used as input for main.jl.
A possible workflow:
1. Define a function q(i,j) on [0,N]^2
2. Write values of q(i,j) into an array and take a real Fourier transform
3. Set the (0,0)th coefficient to be 0 (corresponds to the (1,1)th index)
4. Write Fourier data to output file
"""
function generateIC()
	filename = "initialConditions/vortexPair.csv"


	N = 51 # Require N odd
	q(i,j) = gauss(i,j, [10, 10], 2^2, -6) + gauss(i,j, [65, 40], 10^2, -4) + gauss(i,j, [20, 80], 5^2, 8)

	ICReal = Array{Float64}(undef,(N,N))

	for i=1:N
		for j=1:N
			ICReal[i,j] = q(i,j)
		end
	end


	IC = rfft(ICReal)
	IC[1,1] = 0



	writedlm(filename,  IC, ',')
end

generateIC()

