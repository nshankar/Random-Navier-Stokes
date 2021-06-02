using FFTW
using Random
using Distributions
using Plots
include("helpers.jl")

h = 1
IC, maxFreq = readVelocityIC("IC.csv")
q = curl(IC)
qhat = rfft(q)
triples = computeTriples(maxFreq)
N = 2 * maxFreq + 1

heatmap(q', colors=:grays, aspect_ratio=1, show=true)

for i=1:2000
	t = rand(Gamma(1, h))
	j,k,l = randomTriple(triples)
	println("j: ", j,", k:",  k, ", l:", l, ", tstep:", t)
	evolve!(qhat, j, k, l, t)
	tmp = irfft(qhat, N)
	if mod(i, 10) == 0
		display(heatmap(tmp', colors=:grays, aspect_ratio=1))
		sleep(0.01)
	end
end

