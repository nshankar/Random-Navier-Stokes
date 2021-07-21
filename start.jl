include("main.jl")

function start()
	# Options
	h = 1e-3							# Draws time steps from from Exp(h)		
	cycles = 16							# Number of times to cycle through all of the triples,
										# Each cycle is roughly N^4 ODE calls assuming an N×N grid
	selectTriples = "cyclic"			# One of {"cyclic", "random"}

	fileIC = "initialConditions/gaussianIC.csv"
	fileOutput = "output/long_run10x"			

	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	main(h, Int(cycles), selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
end

start()