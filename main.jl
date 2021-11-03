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

# Helper functions
include("helpers/fileHandling.jl")
include("helpers/mainComputations.jl")
include("helpers/dynamics.jl")
include("helpers/fourierIndexHandling.jl")
include("helpers/animations.jl")

function start()
	# comment out as needed
	Random.seed!(43)


	#Plotting settings
	Plots.scalefontsizes()
	Plots.scalefontsizes(1.2)
	ENV["GKSwstype"]="nul"


	# Options
	h = 5e-5							# Draws time steps from from Exp(h)		
	ncycles = Int(100)					# Number of times to cycle through all of the triples,
										# Each cycle is ~D^2 ODE calls assuming input size D

	fileIC = "initialConditions/smallIC.csv"
	folder = "output/test/"

	#computeVorticityFreqWithViscosity(h, ncycles, fileIC, folder)
	#computeVorticityFreq(h, ncycles, fileIC, folder)
	computeVorticity(folder)

	# Plotting
	animateVorticity(folder, (-5, 5), 20)
	animateLogLogEModes(folder, (-5, 5), 20)
	animateLogEModes(folder, (-5, 5), 20)

	println("All Done!")
end

start()