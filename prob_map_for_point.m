function activation_map = prob_map_for_point(pdf_6D_each, pdf_6D_full, RF_center_points, radii, RF_params, test_point, N_draws)
N_cells = size(RF_center_points,1);

N_grid = 10;
N_pdf = 10;
pdf_x = linspace(0,.2,N_pdf);

activation_map = zeros(N_grid,N_grid,N_draws); %single numbers

for i=1:N_draws
    activation_vec = zeros(N_cells,1);
    activation_ind = zeros(N_cells,1);
    p=zeros(N_grid,N_grid);    
    for c=1:N_cells
        d = pdist2(RF_center_points(c,:),test_point);
        sd = radii(c)/2; %radius is 2 sd of Gauss
        activation_vec(c) = abs((1/(2*pi*sd^2)) * exp(-d^2/(2*sd^2)) + ... %gaussian RF part = mean
            random('Normal',0,RF_params(c).noise_sd)); %noise part = s.d.
        [~, activation_ind(c)] = min(abs(pdf_x - activation_vec(c)));      
    end

    z=1;
    for x=1:N_grid
        for y=1:N_grid
            numerator = prod(pdf_6D_each{z}(activation_ind(1), ...
                activation_ind(2), ...
                activation_ind(3), ...
                activation_ind(4), ...
                activation_ind(5), ...
                activation_ind(6) ...
                ));
            denominator = prod(pdf_6D_full(activation_ind(1), ...
                activation_ind(2), ...
                activation_ind(3), ...
                activation_ind(4), ...
                activation_ind(5), ...
                activation_ind(6) ...
                ));
            activation_map(x,y,i) = numerator/denominator;
            z=z+1;            
        end
    end
end



