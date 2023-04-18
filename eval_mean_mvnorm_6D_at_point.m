function p = eval_mean_mvnorm_6D_at_point(norm_6D, point)
N_grid = size(norm_6D,1);
p = 0;
for x=1:N_grid
    for y=1:N_grid
        p = p + prod(mvnpdf(point, norm_6D(x,y).mu, norm_6D(x,y).sigma));
    end
end
p = p / N_grid^2;