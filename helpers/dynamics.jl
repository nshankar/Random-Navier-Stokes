# This list is not guaranteed to be comprehensive:
# Requires OrdinaryDiffEq
# Requires StaticArrays
# Requires SciMLBase

""" 
    getItpVelocity(qhat)
Use vorticity frequency coefficients qhat to derive a discrete velocity field. 
Then interpolate the velocity values with a quadratic periodic spline.
"""
function getItpVelocity(qhat)
	N = size(qhat)[2]
	vel = biotSavart(qhat)
	grid = LinRange(0, 2*pi, N)
	U = Interpolations.scale(
					interpolate(vel[:,:,1], 
									BSpline(Quadratic(Periodic(OnCell())))),
									grid, grid)
	V = Interpolations.scale(
					interpolate(vel[:,:,2], 
									BSpline(Quadratic(Periodic(OnCell())))),
									grid, grid)
	return U,V
end


"""
    biotSavart(qhat)
Compute the discrete velocity field from vorticity fourier coefficients
using Biot Savart Law, O(N^2 log(N)).
"""
function biotSavart(qhat)
	N = size(qhat)[2]
	maxFreq = Int((N-1)/2)
	velhat = zeros(Complex{Float64}, maxFreq+1, N, 2)
	k = MVector{2,Float64}(0.,0.)
	for k1=0:maxFreq
		for k2=-maxFreq:maxFreq
			k[1] = k1
			k[2] = k2
			qk = qhat[k1+1, mod(k2,N)+1]
			
			if isNonZero(k)
				velhat[k1+1, mod(k2,N)+1,:] = -1im * perp(k)/norm2(k) * qk
			end
		end
	end

	vel = zeros(Float64,N,N,2)
	vel[:,:,1] = irfft(velhat[:,:,1],N)
	vel[:,:,2] = irfft(velhat[:,:,2],N)
	return vel
end


"""
    transport(X, (U,V), t)
Depreciated. Transport a passive scalar by the velocity field (U,V) for time t 
where U,V = getItpVelocity(qhat).
"""
function transport(X, (U,V), t)
	tspan = (0.0, t)
	X_new = zeros(size(X))
	for i=1:size(X)[1]
		prob = ODEProblem(transportODE!, X[i,:], tspan, (U,V))
		sol = solve(prob)
		X_new[i,1] = sol[end][1]
		X_new[i,2] = sol[end][2]
	end
	return rem2pi.(X_new, RoundDown)
end


"""
    transportODE!(dx, x, (U,V), t)
Depreciated. ODE associated with transport(X, (U,V), t).
"""
function transportODE!(dx, x, (U,V), t)
	x[1] = rem2pi(x[1], RoundDown)
	x[2] = rem2pi(x[2], RoundDown)
	dx[1] = U(x[1],x[2])
	dx[2] = V(x[1],x[2])
end


"""
    evolve!(integrator, qhat::Matrix{ComplexF64}, j::SVector{2,Int64}, 
					k::SVector{2,Int64}, l::SVector{2,Int64}, t::Float64, p::MVector{3,Float64})
Take the j,k,l Fourier coefficients and compute their dyanmics for time t. Speed is essential.
"""
function evolve!(integrator, qhat::Matrix{ComplexF64}, j::SVector{2,Int64}, 
					k::SVector{2,Int64}, l::SVector{2,Int64}, t::Float64, p::MVector{3,Float64})
	qj = getFourierCoef(qhat, j)
	qk = getFourierCoef(qhat, k)
	ql = getFourierCoef(qhat, l)
	q0 = SVector{3,ComplexF64}(qj, qk, ql)


	SciMLBase.set_ut!(integrator, q0, 0.)
	add_tstop!(integrator, t)
	integrator.p = p
	solve!(integrator)
	qj, qk, ql = integrator.sol[end]

	setFourierCoef!(qhat, qj, j)
	setFourierCoef!(qhat, qk, k)
	setFourierCoef!(qhat, ql, l)
end


"""
    evolveODE(q0::SVector{3, ComplexF64}, p::MVector{3, Float64}, t::Float64)
ODE associated with evolve!(...). Speed is essential. 
Static vectors are faster than inplace calculation.
"""
function evolveODE(q0::SVector{3, ComplexF64}, p::MVector{3, Float64}, t::Float64)
	C_jl, C_lk, C_kj = p
	qj, qk, ql = q0
	dqj = -conj(C_lk*qk*ql)
	dqk = -conj(C_jl*qj*ql)
	dql = -conj(C_kj*qj*qk)
	return @SVector [dqj, dqk, dql]
end

"""
    getEvolveIntegrator()
Return integrator object initialized to solve evolveODE(...). 
Advantage: the ODE will be initialized infrequently. 
Disadvantage: the integrator appears to cause memory issues.
"""
function getEvolveIntegrator()
	q0 = SVector{3,ComplexF64}(0.,0.,0.)
	tspan = (0., 0.01)
	p = MVector{3,Float64}(0., 0., 0.)
	prob = ODEProblem{false}(evolveODE, q0, tspan, p)
	integrator = init(prob, Tsit5(), save_everystep=false, maxiters = typemax(Int))
	solve!(integrator)
	return integrator
end


"""
    couplingCoef(l, j)
Compute and return C_{l,j}
"""
function couplingCoef(l, j)
	if isNonZero(l) && isNonZero(j)
		return (-j[2]*l[1] + j[1]*l[2])/(4*pi) * (1/norm2(l) - 1/norm2(j))
	end
	return 0
end


### Some functions for specifically length 2 vectors
function norm2(x)
	return x[1]*x[1] + x[2]*x[2]
end

function isNonZero(x)
	if x[1] == 0 && x[2] == 0
		return false
	end
	return true
end

function perp(x)
	x[1], x[2] = -x[2], x[1]
	return x
end


#= moved to fourierIndexHandling.jl
function getFourierCoef(fhat, j)
	N = size(fhat)[2]
	if j[1] >= 0
		fj = fhat[j[1]+1, mod(j[2], N)+1]
	else
		fj = conj(fhat[-j[1]+1, mod(-j[2], N)+1])
	end
	return fj
end

function setFourierCoef!(fhat, fj, j)
	N = size(fhat)[2]
	if j[1] >= 0
		fhat[j[1]+1, mod(j[2], N)+1] = fj
	else
		fhat[-j[1]+1, mod(-j[2], N)+1] = conj(fj)
	end
end
=#