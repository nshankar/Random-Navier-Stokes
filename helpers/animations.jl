# This list is not guaranteed to be comprehensive:
# requires Random
# requires Plots
# requires Printf
# requires StatsPlots
# requires LaTeXStrings


# These functions are good templates, but are not all purpose plotting tools.


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
	mp4(anim, folder*"vorticity.mp4", fps = framerate)
end


# Undergoing rapid changes
function animateLogLogEModes(folder, lims, framerate)
	h, N, ncycles, qData = getVorticityData(folder)
	_, _, _, qhatData = getVorticityFreqData(folder)
	grid = LinRange(0, 2*pi, N)

	qtitlestring = L"\textrm{Vorticity\ Field\ }^{^{}} (h=5 \times 10^{-5})"
	etitlestring = L"\textrm{Flow\ of\ Enstrophy\ &\ Energy}^{^{} }"
	
	plotdims = (500, 500)
	combineddims = (1100, 550)

	anim = @animate for i = 1:ncycles+1
		q = @view qData[:,:,i]
		qhat = @view qhatData[:,:,i]
		p1 = heatmap(grid, grid, q, xticks=(0:pi:2*pi, [L"0" L"\pi" L"2\pi"]), yticks=(0:pi:2*pi, [L"0" L"\pi" L"2\pi"]),
						c = :delta, clims=lims, aspect_ratio=1,
						size=plotdims, title=qtitlestring, colorbar=false)

		logEnergy = computeEVals(qhat, "Energy", true)
		logEnstrophy = computeEVals(qhat, "Enstrophy", true)
		
		p2 = groupedbar([logEnstrophy logEnergy], size=plotdims, title=etitlestring,
					xlabel=L"\log_{3/2}(r)", ylabel=L"\log(1+\overline{E}_r)", xlims=(0,), ylims = (0,18),
					xticks=(0:2:10, [L"0" L"2" L"4" L"6" L"8" L"10"]), yticks=(0:5:15, [L"0" L"5" L"10" L"15"]),
					legend=true, labels=[L"\textrm{Enstrophy}" L"\textrm{Energy}"])
		plot(p1, p2, layout=@layout([A B]), size= combineddims, top_margin=5Plots.mm, bottom_margin=5Plots.mm)
	end
	mp4(anim, folder*"LogLogEModes.mp4", fps = framerate)
end

function animateLogEModes(folder, lims, framerate)
	h, N, ncycles, qData = getVorticityData(folder)
	_, _, _, qhatData = getVorticityFreqData(folder)
	grid = LinRange(0, 2*pi, N)

	qtitlestring = L"\textrm{Vorticity\ Field\ }^{^{}} (h=5 \times 10^{-5})"
	etitlestring = L"\textrm{Flow\ of\ Enstrophy\ &\ Energy}^{^{} }"
	
	plotdims = (500, 500)
	combineddims = (1100, 550)

	anim = @animate for i = 1:ncycles+1
		q = @view qData[:,:,i]
		qhat = @view qhatData[:,:,i]
		p1 = heatmap(grid, grid, q, xticks=(0:pi:2*pi, [L"0" L"\pi" L"2\pi"]), yticks=(0:pi:2*pi, [L"0" L"\pi" L"2\pi"]),
						c = :delta, clims=lims, aspect_ratio=1,
						size=plotdims, title=qtitlestring, colorbar=false)

		logEnergy = computeEVals(qhat, "Energy", false)
		logEnstrophy = computeEVals(qhat, "Enstrophy", false)
		
		p2 = groupedbar([logEnstrophy logEnergy], size=plotdims, title=etitlestring,
					xlabel=L"r", ylabel=L"\log(1+\overline{E}_r)", xlims=(0,70), ylims = (0,18),
					xticks=(0:20:70, [L"0" L"20" L"40" L"60"]), yticks=(0:5:15, [L"0" L"5" L"10" L"15"]),
					legend=true, labels=[L"\textrm{Enstrophy}" L"\textrm{Energy}"])
		plot(p1, p2, layout=@layout([A B]), size= combineddims, top_margin=5Plots.mm, bottom_margin=5Plots.mm)
	end
	mp4(anim, folder*"LogEModes.mp4", fps = framerate)
end


function computeEVals(qhat, E, loglog)
	N = size(qhat)[2]
	maxFreq = Int((N-1)/2)
	if loglog
		scale = 1 + Int(floor(log(1.5, sqrt(2)*maxFreq)))
	else
		scale = Int(floor(sqrt(2)*maxFreq))
	end

	eVals = zeros(scale)
	counts = zeros(scale)
	for i=-maxFreq:maxFreq
		for j=-maxFreq:maxFreq
			if i == j == 0
				continue
			end
			if loglog
				index = 1 + Int(floor(log(1.5, sqrt(i^2+j^2))))
			else
				index = Int(floor(sqrt(i^2+j^2)))
			end

			if E == "Enstrophy"
				eVals[index] = eVals[index] + abs2(getFourierCoef(qhat, [i,j]))
			elseif E == "Energy"
				eVals[index] = eVals[index] + abs2(getFourierCoef(qhat, [i,j]))/(i^2 + j^2)
			else
				return zeros(scale)
			end
			counts[index] = counts[index] + 1
		end
	end


	for i=1:length(counts)
		if counts[i] == 0
			counts[i] = 1
		end
	end
	return log.(ones(scale) + eVals./counts)
end




# This is a temporary solution, but not particularly efficient
function animateVorticityAndEnstrophy(folder, lims, framerate, n_samples)
	h, N, ncycles, qData = getVorticityData(folder)
	_, _, _, qhatData = getVorticityFreqData(folder)
	maxFreq = Int((N-1)/2)

	grid = LinRange(0, 2*pi, N)
	j_indices, k_indices, l_indices = getIndices(maxFreq, n_samples)
	
	# Vectorize this loop
	log_indices = zeros(n_samples)
	for i = 1:n_samples
		log_indices[i] = log(norm2(j_indices[:,i]) + norm2(k_indices[:,i]) + norm2(l_indices[:,i]))
	end

	qtitlestring = Printf.@sprintf "Vorticity Field, h = %.1e, N = %i" h N
	etitlestring = "Log(1 + Enstrophy) Histogram"

	plotdims = (500, 400)
	anim = @animate for i = 1:ncycles+1
		q = @view qData[:,:,i]
		qhat = @view qhatData[:,:,i]
		p1 = heatmap(grid, grid, q,
						c = :delta, clims=lims, 
						size=plotdims, title=qtitlestring)
		
		log_enstrophies = log.(ones(n_samples) + getEnstrophy(qhat, j_indices, k_indices, l_indices))
		#p2 = histogram(log_enstrophies, size=plotdims, title=etitlestring, 
		#				normalize=:pdf, ylims = (0,0.5), xlims = (0, 20), bins=10, legend=false)
		p2 = scatter(log_indices, log_enstrophies, size=plotdims, title="Enstrophy vs Frequency", ylims=(0,15),
						xlabel="log(|j|^2 + |k|^2 + |l|^2)", ylabel="log(1+Enstrophy)", legend=false)
		plot(p1, p2, layout=@layout([A B]), size=(1100, 500), margin=5Plots.mm)
	end
	mp4(anim, folder*"vorticity.mp4", fps = framerate)


end

function getEnstrophy(qhat, j_indices, k_indices, l_indices)
	n_samples = size(j_indices)[2]
	enstrophies = zeros(n_samples)

	for i = 1:n_samples
		qj = getFourierCoef(qhat, j_indices[:,i])
		qk = getFourierCoef(qhat, k_indices[:,i])
		ql = getFourierCoef(qhat, l_indices[:,i])
		enstrophies[i] = abs(qj)^2 + abs(qk)^2 + abs(ql)^2
	end
	return enstrophies
end

function getIndices(maxFreq, n_samples)
	j_indices = zeros(Int64,2,n_samples)
	k_indices = zeros(Int64,2,n_samples)
	l_indices = zeros(Int64,2,n_samples)
	for i=1:n_samples
		j,k,l = getIndex(maxFreq)
		j_indices[:,i] = j
		k_indices[:,i] = k
		l_indices[:,i]= l
	end
	return j_indices, k_indices, l_indices
end

function getIndex(maxFreq)
	j = MVector{2,Int64}(maxFreq+1,maxFreq+1)
	k = MVector{2,Int64}(maxFreq+1,maxFreq+1)
	l = MVector{2,Int64}(maxFreq+1,maxFreq+1)

	# Rejection sampling algorithm: uniformly samples
	# j,k,l in [-maxFreq, maxFreq]^2 such that j + k + l = [0,0]
	for i = 1:2
		while abs(j[i]) > maxFreq || abs(k[i]) > maxFreq || abs(l[i]) > maxFreq
			a = sort(rand(0:3*maxFreq, 2))
			j[i] = a[1] - maxFreq
			k[i] = a[2] - a[1] - maxFreq
			l[i] = 3*maxFreq - a[2] - maxFreq
		end
	end
	return j,k,l
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
	mp4(anim, folder*"all.mp4", fps = framerate)

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

