function cyclicTriple(triples, cycle, ind)
	M = length(triples)
	ind = mod(ind, M^2 - 1) + 1
	res = divrem(cycle[ind], M)

	j = Array{Int64}(undef, 2)
	k = Array{Int64}(undef, 2)
	l = Array{Int64}(undef, 2)


	for i = 1:2
		j[i] = triples[res[i]+1][1]
		k[i] = triples[res[i]+1][2]
		l[i] = triples[res[i]+1][3]
	end
	return j,k,l
end



# Compute a dict of all triples j,k,l in [-N,N] such that j+k+l=0
function computeTriples(N)
	# Consider an immutable dict for speed
	triples = Dict{Int64, Tuple{Int64,Int64,Int64}}()
	counter = 1

	# Consider storing j,k,l in increasing order to save 6x space and some time
	for j = -N:N
		for k = -N:N
			l = - j - k
			if -N <= l <= N
				triples[counter] = (j,k,l)
				counter += 1
			end
		end
	end
	return triples
end


# Uniformly at random selects j,k,l in [-N, N]^2 such that j+k+l=(0,0)
function randomTriple(triples)
	j = Array{Int64}(undef, 2)
	k = Array{Int64}(undef, 2)
	l = Array{Int64}(undef, 2)

	X = rand(1:length(triples), 2)

	for i = 1:2
		j[i] = triples[X[i]][1]
		k[i] = triples[X[i]][2]
		l[i] = triples[X[i]][3]
	end
	return j, k, l
end