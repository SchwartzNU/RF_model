function mv_dist = multi_norm_from_independent_norms(dist_vec)
N = length(dist_vec);

mu = zeros(1,N);
sigma = zeros(N,N);

for i=1:N
    mu(i) = dist_vec{i}.mu;
    sigma(i,i) = dist_vec{i}.sigma;
end

mv_dist.mu = mu;
mv_dist.sigma = sigma;