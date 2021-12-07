"""
    computeTriples(N)
Return dictionary of all integer triples a,b,c in [-N,N] such that a+b+c=0
"""
function computeTriples(N)
	# Consider an array or immutable dict for better space/time efficiency
	triples = Dict{Int64, Tuple{Int64,Int64,Int64}}()
	counter = 1

	# Consider storing j,k,l in increasing order to save 6x space
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


"""
    cyclicTriple(triples, cycle, ind)
For each ind in [1, len(triples)^2], return a unique value of j,k,l in [-N,N]^2 such that j+k+l=(0,0).
The order of the values is determined by cycle, a permutation of length len(triples)^2.
Note: if ind is outside [1, len(triples)^2], it is mapped to its equivalence class in the range.
"""
function cyclicTriple(triples, cycle, ind)
	M = length(triples)
	q_, r_ = divrem(cycle[ind] - 1, M)
	q = q_+1
	r = r_+1

	j = SVector{2,Int64}(triples[q][1],triples[r][1])
	k = SVector{2,Int64}(triples[q][2],triples[r][2])
	l = SVector{2,Int64}(triples[q][3],triples[r][3])

	return j,k,l
end


"""
    randomTriple(triples)
Uniformly at random selects j,k,l in [-N, N]^2 such that j+k+l=(0,0)
"""
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


"""
    getFourierCoef(fhat, j)
Given fhat, a reduced real Fourier transform matrix, 
return the jth Fourier coefficient for any j in [-N,N]^2
"""
function getFourierCoef(fhat, j)
	N = size(fhat)[2]
	if j[1] >= 0
		fj = fhat[j[1]+1, mod(j[2], N)+1]
	else
		fj = conj(fhat[-j[1]+1, mod(-j[2], N)+1])
	end
	return fj
end


"""
    setFourierCoef!(fhat, fj, j)
Given fhat, a reduced real Fourier transform matrix, 
set the jth Fourier coefficient to fj for any j in [-N,N]^2
"""
function setFourierCoef!(fhat, fj, j)
	N = size(fhat)[2]
	if j[1] >= 0
		fhat[j[1]+1, mod(j[2], N)+1] = fj
	else
		fhat[-j[1]+1, mod(-j[2], N)+1] = conj(fj)
	end
end


#=
function myPrint(arr)
	m,n = size(arr)
	for i=1:m
		for j=1:n
			@printf "%.3e+%.3ei " real(arr[i,j]) imag(arr[i,j])
		end
		println()
	end
	println()
end
=#