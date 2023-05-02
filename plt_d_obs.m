close all; clear; clc
load data/d_obs.mat
load data/ObservedData_1Model.mat
load data/results.mat

screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

d_true = [FWD_IP_RC'; FWD_QP_RC'];
d_true = d_true';
d_true = d_true(:);
nlevel = 0.1;
err = nlevel * d_true .* randn(size(d_true));
d_prior = pred_prod(:, :, 1);
d_post = pred_prod(:, :, 5);

noff = 2; nx = 400; ne = 500;
nd = 400;
 
title_names = {'IP-ZZ, offset 1 m', 'IP-ZZ, offset 2 m', ...
               'IP-ZX, offset 1.1 m', 'IP-ZX, offset 2.1 m'};
 
distance = (1:nx)*0.1;
 
figure(1);
set(gcf, 'Position', [0.1*screenWidth 0.1*screenHeight ...
    0.6*screenWidth 0.6*screenHeight])

for i = 1:4
    subplot(2, 2, i)
    start_idx = nd*(i-1)+1;
    end_idx = nd*i;
    idx = start_idx:end_idx;
    hold on
    prior_prct = prctile(d_prior(idx, :)', [5, 95]);
    post_prct = prctile(d_post(idx, :)', [5, 95]);

    patch([distance fliplr(distance)], [prior_prct(1, :) fliplr(prior_prct(2, :))], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])
    patch([distance fliplr(distance)], [post_prct(1, :) fliplr(post_prct(2, :))], [0.3010 0.7450 0.9330], 'EdgeColor', [0.3010 0.7450 0.9330])
    plot(distance, d_true(idx, :), 'color', 'k', 'linewidth', 1.5)
    scatter(distance, d_obs(idx, :), 10, 'r', 'filled');
    if i > 2
        xlabel('Distance (m)');
    else
        xticks([]);
    end
    ylabel('IP (ppm)'); 
    
    title(title_names{i});
    box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title_names = {'QP-ZZ, offset 1 m', 'QP-ZZ, offset 2 m', ...
               'QP-ZX, offset 1.1 m', 'QP-ZX, offset 2.1 m'};
       
figure(2);
set(gcf, 'Position', [0.1*screenWidth 0.1*screenHeight ...
    0.6*screenWidth 0.6*screenHeight])

for i = 1:4
    subplot(2, 2, i)
    start_idx = nd*(i-1)+1+1600;
    end_idx = nd*i+1600;
    idx = start_idx:end_idx;
    hold on
    prior_prct = prctile(d_prior(idx, :)', [5, 95]);
    post_prct = prctile(d_post(idx, :)', [5, 95]);

    patch([distance fliplr(distance)], [prior_prct(1, :) fliplr(prior_prct(2, :))], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])
    patch([distance fliplr(distance)], [post_prct(1, :) fliplr(post_prct(2, :))], [0.3010 0.7450 0.9330], 'EdgeColor', [0.3010 0.7450 0.9330])
    plot(distance, d_true(idx, :), 'color', 'k', 'linewidth', 1.5)
    scatter(distance, d_obs(idx, :), 10, 'r', 'filled');
    if i > 2
        xlabel('Distance (m)');
    else
        xticks([]);
    end
    ylabel('QP (ppm)'); 
    title(title_names{i});
    box on

end