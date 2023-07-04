using JuMP
using Ipopt

model = Model()
"dimension of the problem"
d = 3
"prime modulus"
q = 29
Random.seed!(0)
"b, c: parameters for the optimization problem"
b, c = (rand(0:q-1, d) for _ in 1:2)
@variable(model, a[1:d, 1:d])
# bistochastic constraints:
@constraint(model, a .>= 0)
@constraint(model, a * ones(d) .== ones(d))
# note: writing this as a' * ones(d) does not work:
@constraint(model, ones(d)' * a .== ones(d)')

"""given two elements u, v ∈ ℤ[t]/(t^d+1) represented as vectors,
returns the vector representing the product u⋅v."""
function poly_prod(u, v)
	# coefficient of t^i is ∑ u_j v_{i-j}
	local d = length(u); @assert length(v) == d
	w = zero(u .+ v)
	for i in 1:d, j in 1:d
		k = i+j-1
		x = u[i]*v[j]
		if k > d
			w[k-d] -= x
		else
			w[k] += x
		end
	end
	return w
end

""" poly_prod(u, v, i): same as poly_prod(u, v)[i]. """
function poly_prod(u, v, i)
	s = 0
	local d = length(u); @assert d == length(v)
	for j in 1:d
		k = i-j+1
		if k < 1
			s-= u[j]*v[k+d]
		else
			s+= u[j]*v[k]
		end
	end
	return s
end
"L¹ distance to qℤ^d"
@inline distance_qZ(x::Real) = q*(1/2-abs(mod(x/q,1)-1/2))
@inline distance_qZ(u::Vector) = sum(distance_qZ(u[i]) for i in 1:d)
pp = poly_prod(b, a*c)
dqZ = @NLexpression(model, sum(cos(pp[i]*2π/q) for i in 1:d))


# @NLobjective(model, Min, distance_qZ(poly_prod(b, a * c)))
