# For file IO
using DelimitedFiles

# For main computations
using FFTW
using Random
using Distributions
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using Interpolations

# For animation
using Plots
using StatsPlots
using LaTeXStrings
using ColorSchemes
using Printf

# Helpers
include("helpers/fileHandling.jl")
include("helpers/mainComputations.jl")
include("helpers/dynamics.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/animations.jl")

function main()
	Random.seed!(43)


	# Plotting settings
	Plots.scalefontsizes()
	Plots.scalefontsizes(1.2)
	ENV["GKSwstype"]="nul" # turns off plotting window


	# Options
	h = 5e-5							# Draws time steps from from Exp(h)		
	ncycles = Int(1)					# Number of times to cycle through all of the triples,
											# Each cycle is ~D^2 ODE calls assuming input size D

	fileIC = "initialConditions/vortexPair.csv"
	folder = "output/test/"

	
	lambda(k1, k2) = k1^2 + k2^2
	f(k1, k2) = 1
	
	# Computations
	computeVorticityFreqWithViscosity(h, ncycles, fileIC, folder, lambda, f)
	computeVorticityFreq(h, ncycles, fileIC, folder)
	computeVorticity(folder)

	# Plotting
	animateVorticity(folder, (-5, 5), 20)
	animateLogLogEModes(folder, (-5, 5), 20)
	animateLogEModes(folder, (-5, 5), 20)

	println("All Done!")
end

main()
