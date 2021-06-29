using FFTW
using Random
using Distributions
using OrdinaryDiffEq
using StaticArrays
using DelimitedFiles
using Interpolations
include("helpers/dynamicsHelpers.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/myIO.jl")
include("helpers/myLinAlg.jl")


function main(h, iters, selectTriples, passiveScalars, scalarsCoords, fileIC)
	qhat, maxFreq, vel = readIC(fileIC, "vorticityFreq")

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	if selectTriples == "cyclic"
		cycle = randcycle(length(triples)^2-1)
	end

	# Needs to be fixed
	if passiveScalars == true
		scalarsTraj = zeros(iters, size(scalarsCoords)[1], 2)
		scalarsTraj[1, :, :] = scalarsCoords
		U, V = getItpVelocity(qhat)
	else
		scalarsTraj = nothing
	end

	for i=2:iters
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
		evolve!(qhat, j, k, l, t)

		# Propagate passive scalars
		# Needs to be fixed
		if passiveScalars == true
			U, V = getItpVelocity(qhat)
			scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
		end
	end
	println("All Done!")
end

