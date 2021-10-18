using Plots
using DelimitedFiles
using StaticArrays
include("../../helpers/fileHandling.jl")
include("../../helpers/dynamics.jl")

function energyHistogram(folder)
	h, N, ncycles, data = getVorticityFreqData(folder)
	maxFreq = Int((N-1)/2)

	normalized_energies = Float64[]
	xvals = Float64[]
	qhat = @view data[:,:,ncycles]
	for j1 = -maxFreq:maxFreq
		for j2 = -maxFreq:maxFreq
			for k1 = -maxFreq:maxFreq
				for k2 = -maxFreq:maxFreq
					l1 = j1 + k1
					l2 = j2 + k2
					if abs(l1) < maxFreq && abs(l2) < maxFreq
						j = [j1, j2]
						k = [k1, k2]
						l = [l1, l2]
						C_jl = couplingCoef(j,l)
						C_kl = couplingCoef(l,k)
						C_jk = couplingCoef(k,j)
						qj = getFourierCoef(qhat, j)
						qk = getFourierCoef(qhat, k)
						ql = getFourierCoef(qhat, l)
						energy = (abs(C_kl*qk*ql) + abs(C_jl*qj*ql) + abs(C_jk*qj*ql))/(abs(qj)^2 + abs(qk)^2 + abs(ql)^2)
						push!(xvals, abs(j1)+abs(j2)+ abs(k1)+ abs(k2)+ abs(l1)+ abs(l2))
						push!(normalized_energies, energy)
					end
				end
			end
		end
	end

	histogram(energy)
	png("energies")

	#### Extra!
	# plot random sample of 

end