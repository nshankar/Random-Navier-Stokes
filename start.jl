include("main.jl")

function start()
	# Options
	h = 0.001							# Draws time steps from from Exp(h)		
	cycles = 3							# Number of times to cycle through all of the triples,
										# Each cycle is roughly N^4 ODE calls assuming N×N IC
	eps = 0.01							# Runga Kutta step size (currently not in use)
	selectTriples = "cyclic"			# One of {"cyclic", "random"}

	fileIC = "initialConditions/gaussianIC.csv"
	fileOutput = "output/temp"			

	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	main(h, Int(cycles), eps, selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
end

start()