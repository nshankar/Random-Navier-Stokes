# Workhorse function, takes initial conditions and computes dynamics
function computeVorticityFreq(h, ncycles, fileIC, folder)
	# Read input
	qhat, maxFreq, vel = readIC(fileIC, "vorticityFreq")

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	cycle = randcycle(length(triples)^2)

	# Prepare output file
	mkpath(folder)
	output = open(folder*"vorticityFreq.dat","w")
	writedlm(output, Array{Float64}([h  N ncycles]))
	writedlm(output, qhat, ',')

	# parameters
	p = MVector{3,Float64}(0.,0.,0.)

	println("Computing Vorticity Frequency:")
	for m=1:ncycles
		# Required for memory management
		evolveIntegrator = getEvolveIntegrator()

		for n=1:length(triples)^2
			t = rand(Gamma(1, h))
			# Main dynamics
			j,k,l = cyclicTriple(triples, cycle, n)
			p[1] = couplingCoef(j,l)
			p[2] = couplingCoef(l,k)
			p[3] = couplingCoef(k,j)

			evolve!(evolveIntegrator, qhat, j, k, l, t, p)
		end

		# save output once per cycle
		writedlm(output, qhat, ',')
		println(m) # keep count
		# Enforce garbage collection
		GC.gc()
		sleep(0.001)
	end
	close(output)
end


#Prototype, ideally will be merged with computeVorticityFreq
function computeVorticityFreqWithViscosity(h, ncycles, fileIC, folder, lambda = ((k1,k2)-> k1^2 +k2^2), f = ((k1,k2) -> 1))
	# Read input
	qhat, maxFreq, vel = readIC(fileIC, "vorticityFreq")

	# Prepare random triples
	triples = computeTriples(maxFreq)
	N = 2*maxFreq+1
	cycle = randcycle(length(triples)^2)

	# Prepare output file
	mkpath(folder)
	output = open(folder*"vorticityFreq.dat","w")
	writedlm(output, Array{Float64}([h  N ncycles]))
	writedlm(output, qhat, ',')

	p = MVector{3,Float64}(0.,0.,0.)

	println("Computing Vorticity Frequency:")
	for m=1:ncycles
		# Required for memory management
		evolveIntegrator = getEvolveIntegrator()

		for n=1:length(triples)^2
			t = rand(Gamma(1, h))
			# Main dynamics
			j,k,l = cyclicTriple(triples, cycle, n)
			p[1] = couplingCoef(j,l)
			p[2] = couplingCoef(l,k)
			p[3] = couplingCoef(k,j)
			evolve!(evolveIntegrator, qhat, j, k, l, t, p)
		end

		# Integrate in Viscosity
		k = MVector{2, Int64}(0,0)
		t = rand(Gamma(1,h))
		for k1=0:maxFreq
			for k2=-maxFreq:maxFreq
				k[1] = k1
				k[2] = k2
				qk = getFourierCoef(qhat, k)
				if lambda(k1,k2) == 0
					qk = f(k1,k2)*t + qk
					setFourierCoef!(qhat, qk, k)
				else
					C = qk - f(k1,k2)/lambda(k1,k2)
					qk = C*exp(-lambda(k1,k2)*t) + f(k1,k2)/lambda(k1,k2)
					setFourierCoef!(qhat, qk, k)
				end
			end
		end


		# save output once per cycle
		writedlm(output, qhat, ',')
		println(m) # keep count
		# Enforce garbage collection
		GC.gc()
		sleep(0.001)
	end
	close(output)
end



# Given vorticity frequencies in folder, computes the vorticity
function computeVorticity(folder)
	h, N, ncycles, data = getVorticityFreqData(folder)
	output = open(folder*"vorticity.dat","w")
	writedlm(output, Array{Float64}([h  N ncycles]))

	for i = 1:ncycles+1
		qhat = @view data[:,:,i]
		writedlm(output, irfft(qhat, N), ',')
	end
	close(output)
end

# Given vorticity frequencies in folder, compute the velocity and write to fileoutput
# For now passive scalars are not supported
function computeVelocity(gridSize, passiveScalars, scalarsCoords, folder)
	h, N, ncycles, data = getVorticityFreqData(folder)
	output = open(folder*"velocity.dat","w")
	writedlm(output, Array{Float64}([h  gridSize ncycles]))

	X = repeat(LinRange(0, 2*pi, gridSize), inner=gridSize)
	Y = repeat(LinRange(0, 2*pi, gridSize), outer=gridSize)
	Z = zip(X,Y)

	for i = 1:ncycles+1
		qhat = @view data[:,:,i]
		U, V = getItpVelocity(qhat)
		U_discrete = [U(x,y) for (x,y) in Z]
		V_discrete = [V(x,y) for (x,y) in Z]
		
		writedlm(output, U_discrete, ',')
		writedlm(output, V_discrete, ',')
	end
	close(output)

	# Needs to be fixed
	#=if passiveScalars == true
		scalarsTraj = zeros(iters, size(scalarsCoords)[1], 2)
		scalarsTraj[1, :, :] = scalarsCoords
		U, V = getItpVelocity(qhat)
	else
		scalarsTraj = nothing
	end
	=#

	# Propagate passive scalars
	# Needs to be fixed
	#=if passiveScalars == true
		U, V = getItpVelocity(qhat)
		scalarsTraj[i,:,:] = transport(scalarsTraj[i-1,:,:], (U,V), t)
	end
	=#
end

