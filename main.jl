using FFTW
using Random
using Distributions
using DifferentialEquations
using Plots
using DelimitedFiles
using Interpolations
include("dynamicsHelpers.jl")
include("plottingHelper.jl")
include("fourierIndexHandling.jl")
include("myIO.jl")
include("myLinAlg.jl")


function main(h, iters, eps, selectTriples, passiveScalars, scalarsCoords,
					fileIC, typeIC, viz, limsVorticity, limsVorticityFreq, plotEvery)

	# Error check
	if viz != "velocity" && passiveScalars == true
		println("Visualization type must be velocity to display passive scalars")
		return 1
	end

	qhat, maxFreq, vel = readIC(fileIC, typeIC)

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	if selectTriples == "cyclic"
		cycle = randcycle(length(triples)^2-1)
	end

	# Visualization
	if passiveScalars == true
		scalarsTraj = zeros(iters, size(scalarsCoords)[1], 2)
		scalarsTraj[1, :, :] = scalarsCoords
	else
		scalarsTraj = nothing
	end

	if viz == "velocity" || viz == "all"
		U, V = getItpVelocity(qhat)
	else
		U, V = nothing, nothing
	end
	myplot(1, viz, qhat, (U,V), passiveScalars, scalarsTraj, limsVorticity, limsVorticityFreq)

	for i=2:iters
		t = rand(Gamma(1, h))

		# Main dynamics
		if selectTriples == "cyclic"
			j,k,l = cyclicTriple(triples, cycle, i)
		elseif selectTriples == "random"
			j,k,l = randomTriple(triples)
		else
			return 3
			println("selectTriples ", selectTriples," not recognized.")
		end
		evolve!(qhat, j, k, l, t)

		# Propagate passive scalars
		if passiveScalars == true
			scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
		end

		if mod(i, plotEvery) == 0
		# Visualization 
			if viz == "velocity" || viz == "all"
				U, V = getItpVelocity(qhat)
			end
			myplot(i, viz, qhat, (U,V), passiveScalars, scalarsTraj, limsVorticity, limsVorticityFreq)
		end
	end
	println("All Done!")
end

