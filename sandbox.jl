using FFTW
using Random
using Distributions
using DifferentialEquations
using Plots
using DelimitedFiles
using Interpolations
include("dynamicsHelpers.jl")
include("PlottingHelper.jl")
include("FourierIndexHandling.jl")
include("myIO.jl")
include("myLinAlg.jl")

function bugtest()	
	triples=computeTriples(15)
	println(length(triples)^2)
end


function perms(n, iters)
	for i=1:iters
		v = randperm(n)
	end

end
function myprint(A)
	for i=1:size(A)[1]
		for j=1:size(A)[2]
			print(round(A[i,j], digits=2), " ")
		end
		println()
	end
	println()
end

bugtest()