function activation_vec = activation_for_point(RF_center_points,radii,RF_type,RF_params,p)
%returns a vector with the activation value of each cell - actually it's
%an activation Gaussian for each cell in the non-binary case
%inputs
%RF_center_points: [x, y] for each cell's RF center 
%radii: radius of each RF
%RF_type: 
% 'binary' = 1/area or 0 based on whether point is inside circle of given radius
% 'gauss' = 2D gaussian pdf of the gauss with 2 s.d. = radius
% 'bumpy' = not implemented yet
%RF_params: struct array of parameters for RFs if needed
%  binary: not needed, leave empty
%  gauss: noise_sd = s_d of the additive noise (in units of probabilitiy)
%  bumpy: not implemented yet

N_cells = size(RF_center_points,1);

if strcmp(RF_type,'binary')
    activation_vec = zeros(N_cells,1); %single numbers
else
    activation_vec = cell(N_cells,1); %probability distributions (Gaussian)
end

for i=1:N_cells
    d = pdist2(RF_center_points(i,:),p);
    switch RF_type
        case 'binary'
            if d<radii(i)
                activation_vec(i) = 1/(pi*radii(i)^2);
            else
                activation_vec(i) = 0;
            end
        case 'gauss'
            sd = radii(i)/2; %radius is 2 sd of Gauss
            activation_vec{i} = makedist('Normal', ...
                1/(2*pi*sd^2) * exp(-d^2/(2*sd^2)), ... %gaussian RF part = mean
                RF_params(i).noise_sd);%noise part = s.d.
        case 'bumpy'
    end
end

% %now make map from activation vector
% x_range = [-1, 1];
% y_range = [-1, 1];
% pdf_N = 100; %points in each dimension
% 
% xvals = linspace(x_range(1),x_range(2),pdf_N);
% yvals = linspace(y_range(1),y_range(2),pdf_N);
% [X, Y] = meshgrid(xvals,yvals);
% X_flat = reshape(X,pdf_N^2,1);
% Y_flat = reshape(Y,pdf_N^2,1);
% 
% 
% map_flat = zeros(pdf_N^2,1);
% 
% for i=1:N_cells
%     xvals_shifted = xvals - RF_center_points(i,1);
%     yvals_shifted = yvals - RF_center_points(i,1);
%     
%     curmap = mvnpdf([X_flat, Y_flat], ...
%         repmat(activation_vec{i}.mu,1,2), repmat(activation_vec{i}.sigma,1,2));
%     curmap = reshape(curmap,pdf_N,pdf_N);    
%     imagesc(curmap);
%     colorbar;
%     %pause;
%     keyboard;
% end
% 
% activation_map = reshape(map_flat,pdf_N,pdf_N);
