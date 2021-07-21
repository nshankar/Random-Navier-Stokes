using DelimitedFiles
using Random
using Distributions
using FFTW
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using Interpolations
using LinearAlgebra
include("helpers/dynamics.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/fileHandling.jl")


function main(h, cycles, eps, selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
	# Read input	
	qhat, maxFreq, vel = readIC(fileIC, "vorticityFreq")

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	if selectTriples == "cyclic"
		cycle = randcycle(length(triples)^2)
	end

	# Prepare output file
	output = open(fileOutput,"w")
	writedlm(output, Array{Float64}([h  N]))
	writedlm(output, qhat, ',')

	# Needs to be fixed
	if passiveScalars == true
		scalarsTraj = zeros(iters, size(scalarsCoords)[1], 2)
		scalarsTraj[1, :, :] = scalarsCoords
		U, V = getItpVelocity(qhat)
	else
		scalarsTraj = nothing
	end
	
	evolveIntegrator = getEvolveIntegrator(eps)
	for m=1:cycles
		for n=1:length(triples)^2
			t = rand(Gamma(1, h))

			# Main dynamics
			if selectTriples == "cyclic"
				j,k,l = cyclicTriple(triples, cycle, n)
			elseif selectTriples == "random"
				j,k,l = randomTriple(triples)
			else
				return 1
				println("Triple selection method: ", selectTriples," not recognized.")
			end
			evolve!(evolveIntegrator, qhat, j, k, l, t)

			# save output once per cycle
			if mod(n, N^3) == 0
				writedlm(output, qhat, ',')
			end
		end

		# Propagate passive scalars
		# Needs to be fixed
		if passiveScalars == true
			U, V = getItpVelocity(qhat)
			scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
		end
		# Enforce garbage collection
		GC.gc()
	end
	close(output)
	println("All Done!")
end

