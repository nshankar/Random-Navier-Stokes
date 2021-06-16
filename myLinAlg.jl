# Returns norm squared of a 2d vector
function norm2(x)
	return x[1]^2 + x[2]^2
end

# Returns inner product of 2d vectors
function prod(x,y)
	return x[1]*y[1] + x[2]*y[2]
end

# Returns perp of a 2d vector
function perp(x)
	return [-x[2], x[1]]
end

# Computes curl of 2d vector field using central differences
function curl(u)
	N = size(u)[1]
	delta = 2*pi/N

	u1 = u[:,:,1]
	u2 = u[:,:,2]
	du1_dy = (circshift(u1, (0, 1)) - circshift(u1, (0, -1)))/(2*delta)
	du2_dx = (circshift(u2, (1, 0)) - circshift(u2, (-1, 0)))/(2*delta)
	return du1_dy - du2_dx
end