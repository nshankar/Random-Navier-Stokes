# Requires Plots

function myplot(i, viz, qhat, (U,V), passiveScalars, scalarsTraj, limsVorticity, limsVorticityFreq)
	N = size(qhat)[2]
	plotDims = (450, 500)
	grid = LinRange(0, 2*pi, N)

	if viz == "velocity"
		p = plotVelocity(U,V, plotDims, true)
		if passiveScalars == true
			scalarsCoords = @view scalarsTraj[i,:,:]
			plot!(p, scalarsCoords[:,1], scalarsCoords[:,2], seriestype=:scatter, legend=false )
			#scatter!(scalarsCoords[:,1], scalarsCoords[:,2], legend=false)
			sleep(1/24.)
		end
		display(plot(p))

	elseif viz == "vorticity"
		heatmap(grid, grid, irfft(qhat, N)', colors=:grays, aspect_ratio=1, show=true, clims=limsVorticity,
					size=(600, 600), title="Vorticity Field")
	
	elseif viz == "vorticityFreq"
		heatmap(real.(qhat)', colors=:grays, aspect_ratio=1, show=true, clims=limsVorticityFreq, 
					ylim=(0,N), size=(600, 600), title="Real Part of Fourier Coefs of Vorticity")

	elseif viz == "all"
		p1 = plotVelocity(U,V, plotDims, false)
		p2 = heatmap(grid, grid, irfft(qhat, N)', colors=:grays, aspect_ratio=1, clims=limsVorticity, 
					 size=plotDims, title="Vorticity Field", legend=false)
		p3 = heatmap(real.(qhat)', colors=:grays, clims=limsVorticityFreq, 
						xlim = (1, Int((N+1)/2)), ylim=(1,N), axis=false, size=plotDims, 
						title="Real Part of Fourier Coefs of Vorticity", legend=false)
		display(plot(p1, p2, p3, layout=(1,3), size=(1350, 500)))
	else
		# Save data for later visualization
	end
end


function plotVelocity(U,V, plotDims, disp)
	grid = 10 			# Draws grid^2 vectors
	alpha = 0.5 		# Buffer on x,y limits of velocity field plot
	if disp == true
		plotDims = (600, 600)
	end

	X = repeat(LinRange(0, 2*pi, grid), inner=grid)
	Y = repeat(LinRange(0, 2*pi, grid), outer=grid)
	U_discrete = [U(x,y) for (x,y) in zip(X,Y)]
	V_discrete = [V(x,y) for (x,y) in zip(X,Y)]

	p = quiver(X,Y, quiver=(U_discrete,V_discrete), aspect_ratio=1,
				xlim=(-alpha,2*pi+alpha), ylim=(-alpha,2*pi+alpha),
				size=(600,600), title="Velocity Field")
	return p
end