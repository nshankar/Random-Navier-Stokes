include("main.jl")

function start()
	# Options
	h = 0.1 							# Draws time steps from from Exp(h)				
	iters = 1e3
	eps = 0.001							# Maximum RK4 step size (TBD)	
	selectTriples = "cyclic"			# One of {"cyclic", "random"}

	fileIC = "vorticityFreqIC2.csv" 
	typeIC = "vorticityFreq"			# Initial conditions formatted as one of {"velocity", "vorticityFreq"}

	viz = "vorticity"		   			# Visualize one of {"velocity", "vorticity", "vorticityFreq", "all"}
	limsVorticity = (-3, 3)				# Heatmap limits for vorticity/vorticityFreq visualizations
	limsVorticityFreq = (-100, 100)
	plotEvery = 10						# Plot every N iters

	passiveScalars = false
	scalarsCoords = [1. 3.1; 			
					 2.5 3.;
					 5. 5.]	

	main(h, Int(iters), eps, selectTriples, passiveScalars, scalarsCoords,
			fileIC, typeIC, viz, limsVorticity, limsVorticityFreq, plotEvery)
end

start()