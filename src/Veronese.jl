module Veronese

#
# Utilities.
#
include("util.jl")
export random_coordinates,
    random_vector,
    random_subspace,
    complement_subspace,
    annihilate_subspace,
    coordinate_observation,
    symmetric_product,
    veronese,
    un_veronese,
    subspace_arrangement

#
# Computing veronese rank of subspace arrangements.
#
export subspace_arrangement, incidence_relations, vhs_identifiable

function subspace_arrangement(n, k, d)
    subspaces = zeros(n, d, k)
    tensorized = zeros(n_multipairs(n), n_multipairs(d) * k)

    for i = 1:k
        subspace = random_subspace(n, d)
        subspaces[:, :, i] = subspace
        tensorized[:, (1+n_multipairs(d)*(i-1)):(n_multipairs(d)*i)] = veronese(subspace)
    end

    Vc = annihilate_subspace(tensorized)
    V = annihilate_subspace(Vc)
    return subspaces, V, Vc
end

function incidence_relations(W, Wc, V, Vc)
    n, d = size(V)
    l = size(W, 2)
    nrels = n - d - l + 1

    b = nullspace(Wc' * V, atol = 1e-14)
    a = nullspace(W' * Vc, atol = 1e-14)
    @assert size(b, 2) == 1
    @assert size(a, 2) == nrels

    relations = zeros(d * (n - d), nrels)
    for i = 1:(n-d-l+1)
        relations[:, i] = normalize(reshape(a[:, i] * b', d * (n - d)))
    end
    return relations
end

function generate_summand(subspaces, V, Vc, r)
    n, d, k = size(subspaces)
    sample = subspaces[:, :, rand(1:k)] * random_subspace(d, 1)
    incident = veronese(coordinate_observation(sample, random_coordinates(n, r)))
    W, Wc = complement_subspace(incident)
    relations = incidence_relations(W, Wc, V, Vc)
    relations * relations'
end

function vhs_identifiable(n, k, d, r; verbose = false)
    if r <= d
        return false
    end

    subspaces, V, Vc = subspace_arrangement(n, k, d)

    N, tensorized_dim = size(V)
    dof = tensorized_dim * (N - tensorized_dim)

    dof_per = N - tensorized_dim - binomial(n - r + 2, 2) + 1
    if dof_per <= 0
        return false
    end

    if verbose
        println("truth has dimension $tensorized_dim in $N, $dof DOF")
        println("expected DOF per sample: $dof_per")
        println("expected samples: $(dof / dof_per)")
    end

    var = zeros(dof, dof)

    i = 0
    while i < ceil(dof * 5 / dof_per)
        i += 1
        var += generate_summand(subspaces, V, Vc, r)

        if mod(i, 10) == 0
            c = cond(var)
            if verbose
                println("    $i: $c")
            end
            if c < 1e5
                return true
            end
        end
    end
    return false
end

end
