# For file IO
using DelimitedFiles
include("helpers/fileHandling.jl")

# For main computations
using FFTW
using Random
using Distributions
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using Interpolations
include("main.jl")
include("helpers/dynamics.jl")
include("helpers/fourierIndexHandling.jl")

# For animation
using Plots
using ColorSchemes
using Printf
include("helpers/animations.jl")


function start()
	# Options
	h = 5e-5							# Draws time steps from from Exp(h)		
	ncycles = 10						# Number of times to cycle through all of the triples,
										# Each cycle is roughly N^4 ODE calls assuming an N×N grid

	fileIC = "initialConditions/gaussiansMidRes.csv"
	folder = "output/test/"			

	# passive scalars currently not supported
	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	computeVorticityFreq(h, Int(ncycles), fileIC, folder)
	computeVorticity(folder)
	#computeVelocity(10, passiveScalars, scalarsCoords, folder)


	# Important for plotting behavior
	ENV["GKSwstype"]="nul"
	# Plotting
	animateVorticity(folder, (-5, 5), 12)
	#animateVorticityAndEnstrophy(folder, (-5,5), 6, 1000)
	#animateAll(folder, [(-5, 5), (0, 5)], 10)

	println("All Done!")
end

start()