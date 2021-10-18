# Requires DelimitedFiles

# Read initial conditions
function readIC(fileIC, typeIC)
	if typeIC == "velocity"
		vel, maxFreq = readVelocityIC(fileIC)
		q = curl(vel)
		qhat = rfft(q)
	elseif typeIC == "vorticityFreq"
		qhat, maxFreq = readVorticityFreqIC(fileIC)
		vel = []
	end
	return qhat, maxFreq, vel
end


# Read in velocity initial conditions generated by generateVelocityIC.jl
# Require N odd (for now)
function readVelocityIC(file)
	vel = readdlm(file, ',', Float64)
	N = Int(sqrt(length(vel)/2))
	return reshape(vel, N, N, 2), Int((N - 1)/2)
end

function readVorticityFreqIC(file)
	qhat = readdlm(file, ',', Complex{Float64})
	N = Int( 1/2 * (-1 + sqrt(1 + 8*length(qhat))))
	return qhat, Int((N-1)/2)
end


function getVorticityFreqData(folder)
	f = open(folder*"vorticityFreq.dat", "r")
	data, header = readdlm(f, ',', ComplexF64, header=true)
	close(f)

	h, N, ncycles = split(header[1])
	h = parse(Float64,h)
	N = Int(parse(Float64,N))
	ncycles = Int(parse(Float64, ncycles))

	# This is an unintuitive solution
	dims = (N, Int((N+1)/2), ncycles+1)
	data = reshape(permutedims(data, [2,1]), dims)
	data = permutedims(data, [2,1,3])

	return h, N, ncycles, data
end


function getVorticityData(folder)
	f = open(folder*"vorticity.dat", "r")
	data, header = readdlm(f, ',', Float64, header=true)
	close(f)

	h, N, ncycles = split(header[1])
	h = parse(Float64,h)
	N = Int(parse(Float64,N))
	ncycles = Int(parse(Float64, ncycles))

	dims = (N, N, ncycles+1)
	data = reshape(permutedims(data, [2,1]), dims)
	data = permutedims(data, [2,1,3])

	return h, N, ncycles, data
end



function getVelocityData(folder)
	f = open(folder*"velocity.dat", "r")
	data, header = readdlm(f, ',', Float64, header=true)
	close(f)

	h, gridSize, ncycles = split(header[1])
	h = parse(Float64,h)
	gridSize = Int(parse(Float64,gridSize))
	ncycles = Int(parse(Float64, ncycles))

	dims = (gridSize^2, 2*(ncycles+1))
	data = reshape(data, dims)

	return h, gridSize, ncycles, data
end