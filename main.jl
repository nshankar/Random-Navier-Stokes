using FFTW
using Random
using Distributions
using Plots
include("helpers.jl")

function main(h, iters, file_IC, type_IC, viz, lims_vorticity, lims_vorticityFreq) 
	# Read IC file
	if type_IC == "velocity"
		vel, maxFreq = readVelocityIC(file_IC)
		q = curl(vel)
		qhat = rfft(q)
	elseif type_IC == "vorticityFreq"
		qhat, maxFreq = readVorticityFreqIC(file_IC)
	else
		println("Initial condition type not recognized")
		return 1
	end

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2 * maxFreq + 1

	# Visualization 
	if viz == "velocity" || viz == "all"
		X = repeat(LinRange(0, 2*pi, N), inner=N)
		Y = repeat(LinRange(0, 2*pi, N), outer=N)
		vel = biotSavart(qhat)
		U = reshape(vel[:,:,1], N^2)
		V = reshape(vel[:,:,2], N^2)
	end

	if viz == "velocity"
		quiver(X,Y, quiver=(U,V), arrowscale=0.1, aspect_ratio=1, show=true, 
				xlim=(-0.5,2*pi+0.5), ylim=(-0.5,2*pi+0.5))
	elseif viz == "vorticity"
		heatmap(irfft(qhat, N)', colors=:grays, aspect_ratio=1, show=true, clims=lims_vorticity,
					xlim=(0,N), ylim=(0,N))
	elseif viz == "vorticityFreq"
		heatmap(qhat', colors=:grays, aspect_ratio=1, show=true, clims=lims_vorticityFreq, ylim=(0,N))
	elseif viz == "all"
		p1 = quiver(X,Y, quiver=(U,V), arrowscale=0.1, aspect_ratio=1, 
						xlim=(-0.5,2*pi+0.5), ylim=(-0.5,2*pi+0.5))
		p2 = heatmap(irfft(qhat, N)', colors=:grays, aspect_ratio=1, clims=lims, xlim=(0,N), ylim=(0,N))
		p3 = heatmap(qhat', colors=:grays, aspect_ratio=1, clims=lims, ylim=(0,N))
		display(plot(p1, p2, p3, layout=(1,3)))#, size=(1000, 400)))
	else
		println("Visualization type not recognized")
		return 2
	end
	sleep(2)

	for i=1:iters
		t = rand(Gamma(1, h))
		j,k,l = randomTriple(triples)
		# Verbose
		# println("j: ", j,", k:",  k, ", l:", l, ", tstep: ", t)
		evolve!(qhat, j, k, l, t)
		
		# Update plot every 10 iters
		if mod(i, 10) == 0
			# Visualization 
			if viz == "velocity" || viz == "all"
				vel = biotSavart(qhat)
				U = reshape(vel[:,:,1], N^2)
				V = reshape(vel[:,:,2], N^2)
			end

			if viz == "velocity"
				quiver(X,Y, quiver=(U,V), arrowscale=0.1, aspect_ratio=1, show=true, xlim=(-0.5,2*pi+0.5), 
						ylim=(-0.5,2*pi+0.5))
			elseif viz == "vorticity"
				heatmap(irfft(qhat, N)', colors=:grays, aspect_ratio=1, show=true, clims=lims_vorticity,
						xlim=(0,N), ylim=(0,N))
			elseif viz == "vorticityFreq"
				heatmap(qhat', colors=:grays, aspect_ratio=1, show=true, clims=lims_vorticityFreq,
						ylim=(0,N))
			else # viz == "all"
				p1 = quiver(X,Y, quiver=(U,V), arrowscale=0.1, aspect_ratio=1, 
							xlim=(-0.5,2*pi+0.5), ylim=(-0.5,2*pi+0.5))
				p2 = heatmap(irfft(qhat, N)', colors=:grays, aspect_ratio=1, clims=lims_vorticity, 
								ylim=(-1,N+1))
				p3 = heatmap(qhat', colors=:grays, aspect_ratio=1, clims=lims_vorticityFreq, ylim=(0,N))
				display(plot(p1, p2, p3, layout=(1,3)))#, size=(1000, 400)))
			end
			println(".")
		end


	end
end

