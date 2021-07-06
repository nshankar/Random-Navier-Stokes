include("main.jl")

function start()
	# Options
	h = 0.01 							# Draws time steps from from Exp(h)		
	iters = 100*157609					# Enough for a 1000*h length simulation
	selectTriples = "cyclic"			# One of {"cyclic", "random"}

	fileIC = "initialConditions/vorticityFreqIC3.csv"
	fileOutput = "output/output3"		# terrible naming convention, fix me

	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	main(h, Int(iters), selectTriples, passiveScalars, scalarsCoords, fileIC, fileOutput)
end

start()