using DelimitedFiles
using Random
using Distributions
using FFTW
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using Interpolations
include("helpers/dynamicsHelpers.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/myIO.jl")
include("helpers/myLinAlg.jl")


function main(h, iters, selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
	# Read input	
	qhat, maxFreq, vel = readIC(fileIC, "vorticityFreq")


	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	if selectTriples == "cyclic"
		cycle = randcycle(length(triples)^2 - 1) # note: why -1? // Remark: there is an issue that needs to be fixed here
	end

	# Prepare output file
	output = open(fileOutput,"w")
	writedlm(output, Array{Float64}([h  N]))


	# Needs to be fixed
	if passiveScalars == true
		scalarsTraj = zeros(iters, size(scalarsCoords)[1], 2)
		scalarsTraj[1, :, :] = scalarsCoords
		U, V = getItpVelocity(qhat)
	else
		scalarsTraj = nothing
	end


	
	evolveIntegrator = getEvolveIntegrator()
	for i=1:iters
		t = rand(Gamma(1, h))

		# Main dynamics
		if selectTriples == "cyclic"
			j,k,l = cyclicTriple(triples, cycle, i)
		elseif selectTriples == "random"
			j,k,l = randomTriple(triples)
		else
			return 1
			println("Triple selection method: ", selectTriples," not recognized.")
		end
		evolve!(evolveIntegrator, qhat, j, k, l, t)

		# Propagate passive scalars
		# Needs to be fixed
		if passiveScalars == true
			U, V = getItpVelocity(qhat)
			scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
		end

		# save output once per cycle
		if mod(i, 157609) == 1
			writedlm(output, qhat, ',')
		end
	end
	close(output)
	println("All Done!")
end

