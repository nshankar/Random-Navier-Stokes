using Plots
using DelimitedFiles
using StaticArrays
include("../../helpers/dynamics.jl")
include("../../helpers/fourierIndexHandling.jl")


# for each k, find the max C_jl such that j + l = k
function computeCouplingCoefsSupremumStandard(N)
	maxFreq = Int((N-1)/2)
	l1 = 0
	l2 = 0
	A = zeros(N,N)
	for k1 = -maxFreq:maxFreq
		for k2 = -maxFreq:maxFreq
			max = 0
			for j1 = -maxFreq:maxFreq
				for j2 = -maxFreq:maxFreq
					l1 = k1 - j1
					l2 = k2 - j2
					if abs(l1) < maxFreq && abs(l2) < maxFreq
						if abs(couplingCoef([j1, j2], [l1, l2])) > max
							max = abs(couplingCoef([j1, j2], [l1, l2]))
						end
					end
				end
			end

			A[maxFreq+k1+1, maxFreq+k2+1] = max
		end
	end
	return heatmap(A, c=:heat, axis=false, title="N = "*string(N), xticks=false, yticks=false)
	#png("couplingCoefsSupremumStandard")
end

# for each j find the max value of C_{jk}
function computeCouplingCoefsSupremum(N)
	maxFreq = Int((N-1)/2)
	j = MVector{2, Int64}(0,0)
	k = MVector{2, Int64}(0,0)
	A = zeros(N,N)
	for j1 = -maxFreq:maxFreq
		j[1] = j1
		for j2 = -maxFreq:maxFreq
			j[2] = j2
			max = 0
			for k1 = -maxFreq:maxFreq
				k[1] = k1
				for k2 = -maxFreq:maxFreq
					k[2] = k2

					if abs(couplingCoef(j,k)) > max
						max = abs(couplingCoef(j,k))
					end

				end
			end

			A[maxFreq+j[1]+1, maxFreq+j[2]+1] = max
		end
	end

	return heatmap(A, c=:delta, axis=false, title="N = "*string(N), xticks=false, yticks=false)		
end


function plotMany()
	Ns = [11, 51, 101, 151, 201, 251]

	p1 = computeCouplingCoefsSupremumStandard(Ns[1])
	p2 = computeCouplingCoefsSupremumStandard(Ns[2])
	p3 = computeCouplingCoefsSupremumStandard(Ns[3])
	p4 = computeCouplingCoefsSupremumStandard(Ns[4])
	p5 = computeCouplingCoefsSupremumStandard(Ns[5])
	p6 = computeCouplingCoefsSupremumStandard(Ns[6])

	title = plot(title="Coupling Coefs Max Values",framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0Plots.px)


	plot(title, p1, p2, p3, p4, p5, p6, layout=@layout([A{0.05h}; [B C D]; [E F G]]), size=(1350, 550))
	png("couplingCoefsSizeStandardComparisons")
end



function computeGrowthRate()
	Ns = [11, 101, 201, 301, 401, 501, 601, 701, 801, 901, 1001]
	data = zeros(length(Ns))
	# uses shortcut that j = [1,2] is essentially maximal
	for i = 1:length(Ns)
		N = Ns[i]
		maxFreq = Int((N-1)/2)
		j = SVector{2, Int64}(1,2)
		k = MVector{2, Int64}(0,0)
		max = 0
		for k1 = -maxFreq:maxFreq
			k[1] = k1
			for k2 = -maxFreq:maxFreq
				k[2] = k2

				if abs(couplingCoef(j,k)) > max
					max = abs(couplingCoef(j,k))
				end
			end
		end
		data[i] = max
	end
	scatter(Ns, data, xlabel="N", ylabel="Max Coef", legend = false, title="Coupling Coefficient Growth Rate")
	png("growthRate")
	#output = open("extras/couplingCoefsAnalysis/growthRate.dat","w")
	#writedlm(output, Ns, ',')
	#writedlm(output, data, ',')
	#close(output)
end