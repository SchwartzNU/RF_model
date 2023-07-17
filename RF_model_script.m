noise_levels = [1E-2, 5E-3, 1E-3, 5E-4, 1E-4, 5E-5, 1E-5];
radius_levels = 1:.25:3;

N_noise = length(noise_levels);
N_radius = length(radius_levels);

for n=1:N_noise
    for r = 1:N_radius        
        tic;
        fprintf('Noise = %d, Radius = %d\n', noise_levels(n), radius_levels(r));
        fname = sprintf('posterior_roise_%d_radius_%d.mat', noise_levels(n), radius_levels(r));
        [RGC_X,RGC_Y] = Hex_grid_generator(3,3);
        RGC_X = RGC_X(1:end-1)-1;
        RGC_Y = RGC_Y(1:end-1)-1;

        radius = radius_levels(r);

        N_cells = length(RGC_X);
        RF_center_points = [RGC_X, RGC_Y];
        radii = ones(N_cells,1)*radius;

        param = struct;
        param.noise_sd = noise_levels(n);
        RF_params = repmat(param,N_cells,1);

        posterior = posterior_map_for_point(RF_center_points, radii, RF_params, [-1 1], [-1 1], [.25 .5], 60, 10);

        save(fname, 'posterior');
        elapsed = toc;
        fprintf('%d seconds \n', elapsed);
    end
end

%%
point_location = [.25 .5];
N_grid = 50;
rangeX = [-1 1];
rangeY = [-1 1];
xvals = linspace(rangeX(1),rangeX(2),N_grid);
yvals = linspace(rangeY(1),rangeY(2),N_grid);

posterior_gauss_fit = cell(N_noise, N_radius);
posterior_gauss_sigX = zeros(N_noise, N_radius);
posterior_gauss_sigY = zeros(N_noise, N_radius);
posterior_gauss_muX = zeros(N_noise, N_radius);
posterior_gauss_muY = zeros(N_noise, N_radius);
posterior_gauss_ang = zeros(N_noise, N_radius);
z=1;
for n=1:N_noise
    coverage_factor = zeros(N_radius,N_grid^2);
    for r = 1:N_radius        
        fprintf('Noise = %d, Radius = %d\n', noise_levels(n), radius_levels(r));
        fname = sprintf('posterior_roise_%d_radius_%d.mat', noise_levels(n), radius_levels(r));
        [RGC_X,RGC_Y] = Hex_grid_generator(3,3);
        RGC_X = RGC_X(1:end-1)-1;
        RGC_Y = RGC_Y(1:end-1)-1;

        z = 1;
%         for x=1:N_grid
%             for y=1:N_grid
%                 for i=1:length(RGC_X)            
%                     d = pdist2([xvals(x), yvals(y)], [RGC_X(i), RGC_Y(i)]);
%                     if d <= radius_levels(r)                        
%                         coverage_factor(r,z) = coverage_factor(r,z) + 1;                        
%                     end
%                 end
%                 z=z+1;
%             end
%         end
%        fprintf('Radius = %d, Coverage factor = %d \n', radius_levels(r), mean(coverage_factor(r,:)));
        load(fname);
        figure(z);
        subplot(1,2,1);
        imagesc(posterior);
        [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(xvals,yvals,posterior);
        posterior_gauss_fit{n,r} = zfit;
        posterior_gauss_sigX(n,r) = fitresult(3);
        posterior_gauss_sigY(n,r) = fitresult(4);
        posterior_gauss_muX(n,r) = fitresult(5);
        posterior_gauss_muY(n,r) = fitresult(6);
        posterior_gauss_ang(n,r) = fitresult(2);
        subplot(1,2,2);
        imagesc(zfit);
    end 
end

%% compute coverage 
N_grid = 50;
rangeX = [-1 1];
rangeY = [-1 1];
xvals = linspace(rangeX(1),rangeX(2),N_grid);
yvals = linspace(rangeY(1),rangeY(2),N_grid);

for n=1:N_noise
    coverage_factor = zeros(N_radius,N_grid^2);
    for r = 1:N_radius 
        z=1;
        for x=1:N_grid
            for y=1:N_grid
                for i=1:length(RGC_X)            
                    d = pdist2([xvals(x), yvals(y)], [RGC_X(i), RGC_Y(i)]);
                    if d <= radius_levels(r)                        
                        coverage_factor(r,z) = coverage_factor(r,z) + 1;                        
                    end
                end
                z=z+1;
            end
        end
    end
end

%% make RGC grid
[RGC_X,RGC_Y] = Hex_grid_generator(3,3);
RGC_X = RGC_X(1:end-1)-1;
RGC_Y = RGC_Y(1:end-1)-1;
%6 RFS in a perfect hex grid centered at 0,0
%Distance between centers for neighbors is 2 units

radius = 2;

N_cells = length(RGC_X);
RF_center_points = [RGC_X, RGC_Y];
radii = ones(N_cells,1)*radius;

param = struct;
params.noise_sd = 1E-5;
RF_params = repmat(params,N_cells,1);

%posterior_noise_5 = posterior_map_for_point(RF_center_points, radii, RF_params, [-1 1], [-1 1], [.25 .5], 40, 10)

%%
for r = 1:length(radius_levels)
    figure(r);
    % plot the centers
    scatter(RGC_X, RGC_Y, 50, 'k+');
    ax(r) = gca;
    xlim(ax(r), [-4 4]);
    ylim(ax(r), [-4 4]);

    % make circle rois
    colors = 'rgbcmy';
    for i=1:6
        RF(i) = drawcircle(ax(r),'center',[RGC_X(i),RGC_Y(i)],...
            'color',colors(i),...
            'FaceAlpha', 0.2, ...
            'InteractionsAllowed','none',...
            'linewidth',1,...
            'radius',radius_levels(r));
    end
    ch = get(ax(r), 'children');   
    circles = [];
    for c=1:length(ch)
        if strcmp(class(ch(c)), 'images.roi.Circle')
            circles = [circles, get(ch(c), 'vertices')];
        end
    end
    circles = circles .* factor(r);
    fname = sprintf('Circles_radius_%d.txt', radius_levels(r));
    dlmwrite(fname, circles);
end

%% make RF examples with different SNR

f = fspecial('gaussian',[600, 600], 100);
f_noise = zeros(600, 600, N_noise);
for n=1:N_noise
    f_noise(:,:,n) = f + random('Normal',0,noise_levels(n)*noise_scale,[600 600]);
    %imagesc(f_noise(:,:,n));
    fname = sprintf('RF_noise_%d.txt', n);
    dlmwrite(fname, f_noise(:,:,n));
end

