function animateVorticity(folder, lims, framerate)
	h, N, ncycles, data = getVorticityData(folder)
	grid = LinRange(0, 2*pi, N)
	titlestring = Printf.@sprintf "Vorticity Field, h = %.1e, N = %i" h N

	anim = @animate for i = 1:ncycles+1
		q = @view data[:,:,i]
		heatmap(grid, grid, q,
				c = :delta,	aspect_ratio=1, clims=lims, 
				show=false, title=titlestring)
	end
	gif(anim, folder*"vorticity.gif", fps = framerate)
end


function animateAll(folder, lims, framerate)
	h, N, ncycles, qSeries = getVorticityData(folder)
	h, N, ncycles, qhatSeries = getVorticityFreqData(folder)
	h, gridSize, ncycles, velSeries = getVelocityData(folder)

	grid = LinRange(0, 2*pi, N)
	X = repeat(LinRange(0, 2*pi, gridSize), inner=gridSize)
	Y = repeat(LinRange(0, 2*pi, gridSize), outer=gridSize)
	# buffer at edge of plot
	alpha = 0.5 
	titlestring = Printf.@sprintf "h = %.1e, N = %i" h N
	plotdims = (450, 500)

	anim = @animate for i = 1:ncycles+1
		q = @view qSeries[:,:,i]
		qhat = @view qhatSeries[:,:,i]
		U_discrete = velSeries[:, 2*i-1]
		V_discrete = velSeries[:, 2*i]

		# Vorticity Field
		p1 = heatmap(grid, grid, q,
				c = :delta,	aspect_ratio=1, clims=lims[1], 
				show=false, title="Vorticity", legend=false)

		# Vorticity Fourier
		p2 = heatmap(log.(abs.(qhat) + ones(size(qhat))), c=:blues, clims=lims[2], 
						xlim = (1, N), ylim=(1, Int((N+1)/2)), axis=false, 
						size=plotdims, 	title="Log Vorticity Freq", legend=false)
		
		# Velocity Field
		# What's with the negated transpose?
		p3 = quiver(Y, X, quiver=(-V_discrete,-U_discrete), aspect_ratio=1,
				xlim=(-alpha,2*pi+alpha), ylim=(-alpha,2*pi+alpha),
				size=plotdims, title="Velocity Field")

		# Title
		title = plot(title=titlestring,framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0Plots.px)

		plot(title, p1, p2, p3, layout=@layout([A{0.05h}; [B C D]]), size=(1350, 550))
	end
	gif(anim, folder*"all.gif", fps = framerate)

end

function plotVelocity(U,V, plotDims)
	grid = 10 			# Draws grid^2 vectors
	alpha = 0.5 		# Buffer on x,y limits of velocity field plot

	X = repeat(LinRange(0, 2*pi, grid), inner=grid)
	Y = repeat(LinRange(0, 2*pi, grid), outer=grid)
	U_discrete = [U(x,y) for (x,y) in zip(X,Y)]
	V_discrete = [V(x,y) for (x,y) in zip(X,Y)]

	p = quiver(X,Y, quiver=(U_discrete,V_discrete), aspect_ratio=1,
				xlim=(-alpha,2*pi+alpha), ylim=(-alpha,2*pi+alpha),
				size=plotDims, title="Velocity Field")
	return p
end

