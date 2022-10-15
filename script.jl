using Veronese

# # adjust the parameters n, k, d here
# arrangement_well_posed(21, 7, 3, 4)

n, k, d, r = 10, 3, 3, 5

data, V, Vc = subspace_arrangement(n, k, d)
N, tensorized_dim = size(V)
dof = tensorized_dim * (N - tensorized_dim)
var = zeros(dof, dof)
sample = data[:, :, rand(1:k)] * random_subspace(d, 1)
incident = veronese(coordinate_observation(sample, random_coordinates(n, r)))
W, Wc = complement_subspace(incident)
