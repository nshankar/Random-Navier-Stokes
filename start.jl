include("main.jl")

function start()
	# Options
	h = 0.1 							# Draws time steps from from Exp(h)			
	iters = 1e4
	selectTriples = "cyclic"			# One of {"cyclic", "random"}

	fileIC = "initialConditions/vorticityFreqIC2.csv"

	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	main(h, Int(iters), selectTriples, passiveScalars, scalarsCoords,fileIC)
end

start()