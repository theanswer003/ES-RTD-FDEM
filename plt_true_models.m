clear; close all; clc
addpath(genpath('./utils'));
addpath(genpath('./forward'));

load data/EC_TrueModel.mat
load data/MS_TrueModel.mat
load data/synthetic_models.mat
load data/results.mat


screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

nx = 400; nz = 40;
EC_true = squeeze(EC_TrueModel);
MS_true = squeeze(MS_TrueModel);

post_EC_mu = mean(post_EC, 3);
post_MS_mu = mean(post_MS, 3);
post_EC_std =std(post_EC, 0, 3);
post_MS_std = std(post_MS, 0, 3);

prior_EC_mu = mean(EC_prior, 3);
prior_MS_mu = mean(MS_prior, 3);

minEC = min(EC_true(:));
maxEC = max(EC_true(:));
minMS = min(MS_true(:));
maxMS = max(MS_true(:));

figure(1);
set(gcf, 'Position', [0.1*screenWidth 0.1*screenHeight ...
    0.8*screenWidth 0.4*screenHeight])
subplot(221)
imagesc([0, nx]*0.1, [0, nz]*0.1, EC_true'); 
xticks([]); ylabel('Depth (m)')
title('True EC')
caxis([minEC, maxEC]); c1 = colorbar; ylabel(c1, 'S/m');
dim = [.09+0.01 .87 .1 .1];
h = annotation('textbox', dim, 'String', '(a)', ...
    'fontsize', 14, 'fontweight', 'bold','FitBoxToText','on');
h.LineStyle = 'none';

subplot(222)
imagesc([0, nx]*0.1, [0, nz]*0.1, MS_true'); 
xticks([]); ylabel('Depth (m)')
title('True MS')
caxis([minMS, maxMS]); c2 = colorbar; ylabel(c2, 'SI');
dim = [.53+0.01 .87 .1 .1];
h = annotation('textbox', dim, 'String', '(b)', ...
    'fontsize', 14, 'fontweight', 'bold','FitBoxToText','on');
h.LineStyle = 'none';

subplot(223)
imagesc([0, nx]*0.1, [0, nz]*0.1, prior_EC_mu'); 
xlabel('Distance (m)'); ylabel('Depth (m)'); xticks(5:5:35);
title('Prior mean of EC')
caxis([minEC, maxEC]); c3 = colorbar ; ylabel(c3, 'S/m');
dim = [.09+0.01 .4 .1 .1];
h = annotation('textbox', dim, 'String', '(c)', ...
    'fontsize', 14, 'fontweight', 'bold','FitBoxToText','on');
h.LineStyle = 'none';

subplot(224)
imagesc([0, nx]*0.1, [0, nz]*0.1, prior_MS_mu'); 
xlabel('Distance (m)'); ylabel('Depth (m)'); xticks(5:5:35);
title('Prior mean of MC')
caxis([minMS, maxMS]); c4 = colorbar ; ylabel(c4, 'SI');
dim = [.53+0.01 .4 .1 .1];
h = annotation('textbox', dim, 'String', '(d)', ...
    'fontsize', 14, 'fontweight', 'bold','FitBoxToText','on');
h.LineStyle = 'none';
colormap('jet')

