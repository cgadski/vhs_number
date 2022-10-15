using LinearAlgebra
using Random: randperm

random_coordinates(n, k) = BitVector(
	sum(Matrix(I, n, n)[:, randperm(n)[1:k]], dims = 2)[:, 1]
	)

function random_vector(n::Int) 
	v = randn(n)
	(1 / norm(v)) * v
end

random_subspace(n::Int, d::Int) = qr(reshape(randn(n * d), n, d)).Q[:, 1:d]

function complement_subspace(m)
	fact = qr(m)
	rank = sum(abs.(diag(fact.R)) .> 1e-15)
	fact.Q[:, 1:rank], fact.Q[:, (rank + 1):end]
end

annihilate_subspace(m) = nullspace(transpose(m))

"Incident subspace associated with a coordinate observation."
function coordinate_observation(vector, observed)
	n = length(vector)
	k = n - sum(observed)

	incident::Matrix{Float64} = zeros((length(vector), k + 1))
	begin
		idx = 1
		for i = 1:n
			if !observed[i]
				incident[i, idx] = 1
				idx += 1
			end
		end

		incident[:, end] = normalize([!observed[i] ? 0 : vector[i] for i in 1:n])
	end

	return incident
end

n_multipairs(n::Integer) = binomial(n + 1, 2)

function symmetric_product(v::Vector{Float64}, w::Vector{Float64})
	n = length(v)
	o = Vector{Float64}(undef, n_multipairs(n))
	k = 1
	for i in 1:n
		for j in i:n
			o[k] = i == j ? v[i] * w[i] : (v[i] * w[j] + v[j] * w[i])
			k += 1
		end
	end
	return o
end

veronese(v::Vector{Float64}) = symmetric_product(v, v)

"The symmetric tensor power for vectors or linear maps."
function veronese(m::Matrix{Float64})
	n, d = size(m)
	columns::Matrix{Float64} = zeros((n_multipairs(n), n_multipairs(d)))
	k = 1
	for i in 1:d
		for j in i:d
			columns[:, k] = symmetric_product(m[:, i], m[:, j])
			k += 1
		end
	end
	return columns
end

function vector_to_symmetric(v::Vector{Float64}, n)
	m::Matrix{Float64} = Matrix(undef, n, n)
	k = 1
	for i in 1:n
		for j in i:n
			if i == j
				m[i, j] = v[k]
			else
				m[i, j] = m[j, i] = 0.5 * v[k]
			end
			k += 1
		end
	end
	return Symmetric(m)
end

function un_veronese(v::Vector{Float64}, n::Int)
	mat = vector_to_symmetric(v, n)
	lambda, v = eigen(mat, n:n)
	return sqrt(lambda[1]) * v[:, 1]
end

