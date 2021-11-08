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
	ENV["GKSwstype"]="nul" # turns off plotting window


	# Options
	h = 5e-5							# Draws time steps from from Exp(h)		
	ncycles = Int(5)					# Number of times to cycle through all of the triples,
										# Each cycle is ~D^2 ODE calls assuming input size D

	fileIC = "initialConditions/gaussianIC.csv"
	folder = "output/test/"

	N = 101
	lambda(k1, k2) = Int((k1^2 + k2^2)*((abs(k1) <= 4 && abs(k2) <= 4) || (abs(k1) >= N//2-4 || abs(k2)>= N//2 -4)))
	f(k1, k2) = Int(abs(k1) <= 4 && abs(k2) <= 4)
	computeVorticityFreqWithViscosity(h, ncycles, fileIC, folder, lambda, f)
	#computeVorticityFreq(h, ncycles, fileIC, folder)
	computeVorticity(folder)

	# Plotting
	animateVorticity(folder, (-5, 5), 20)
	animateLogLogEModes(folder, (-5, 5), 20)
	animateLogEModes(folder, (-5, 5), 20)

	println("All Done!")
end

start()