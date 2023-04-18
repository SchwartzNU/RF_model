function [pdf_6D_full, activation_map, entropy_map, max_prob_map] = mv_gauss_map(RF_center_points, radii, RF_params, rangeX, rangeY, N)

N_pdf = 10;
pdf_x = linspace(0,.2,N_pdf);
[g1, g2, g3, g4, g5, g6] = ndgrid(pdf_x);

g1_flat = reshape(g1,N_pdf^6,1);
g2_flat = reshape(g2,N_pdf^6,1);
g3_flat = reshape(g3,N_pdf^6,1);
g4_flat = reshape(g4,N_pdf^6,1);
g5_flat = reshape(g5,N_pdf^6,1);
g6_flat = reshape(g6,N_pdf^6,1);

%map = zeros(N,N);
xvals = linspace(rangeX(1),rangeX(2),N);
yvals = linspace(rangeY(1),rangeY(2),N);
[X, Y] = meshgrid(xvals,yvals);
entropy_map = zeros(N,N);
activation_map = zeros(N,N,6);
max_prob_map = zeros(N,N);
pdf_6D_each = cell(N^2,1);
posterior = zeros(N,N);
z = 1;
for i=1:N
    i
    for j=1:N
        activation_vec = activation_for_point(RF_center_points,radii,'gauss',RF_params,[X(i,j), Y(i,j)]);
        mv_dist = multi_norm_from_independent_norms(activation_vec);
        
        pdf_flat = mvnpdf([g1_flat, g2_flat, g3_flat, g4_flat, g5_flat, g6_flat],...
            mv_dist.mu, mv_dist.sigma);

        pdf_6D = reshape(pdf_flat,[N_pdf,N_pdf,N_pdf,N_pdf,N_pdf,N_pdf]);
        pdf_6D_each{z} = pdf_6D ./ sum(pdf_6D(:)); %normalize'

        [max_val, ind] = max(pdf_6D(:));
        max_prob_map(i,j) = max_val;        
        [p1, p2, p3, p4, p5, p6] = ind2sub(size(pdf_6D),ind);
        activation_map(i,j,1) = pdf_x(p1);
        activation_map(i,j,2) = pdf_x(p2);
        activation_map(i,j,3) = pdf_x(p3);
        activation_map(i,j,4) = pdf_x(p4);
        activation_map(i,j,5) = pdf_x(p5);
        activation_map(i,j,6) = pdf_x(p6);
        entropy_map(i,j) = entropy_from_pdf(pdf_6D);
        
        if i==1 && j==1
            pdf_6D_full = pdf_6D;
        else
            pdf_6D_full = pdf_6D_full + pdf_6D;
        end
        z=z+1;
    end
end
pdf_6D_full = pdf_6D_full./N^2;

% z = 1;
% for i=1:N
%     i
%     for j=1:N
%         posterior(i,j) = nansum(pdf_6D_each{z}./pdf_6D_full,'all');
%         z=z+1;
%     end
% end
