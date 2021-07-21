using DelimitedFiles
using Random
using Distributions
using FFTW
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using Interpolations
include("helpers/dynamics.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/fileHandling.jl")


function main(h, cycles, selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
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
	
	p = MVector{3,Float64}(0.,0.,0.)
	for m=1:cycles
		# attempt to solve memory issues
		evolveIntegrator = getEvolveIntegrator()

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
			p[1] = couplingCoef(j,l)
			p[2] = couplingCoef(l,k)
			p[3] = couplingCoef(k,j)
			evolve!(evolveIntegrator, qhat, j, k, l, t, p)

		end

		# Propagate passive scalars
		# Needs to be fixed
		if passiveScalars == true
			U, V = getItpVelocity(qhat)
			scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
		end

		# save output once in a while
		writedlm(output, qhat, ',')
		# Enforce garbage collection
		GC.gc()
		sleep(0.001)
	end
	close(output)
	println("All Done!")
end

