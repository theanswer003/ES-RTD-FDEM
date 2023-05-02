clear; close all; clc
addpath(genpath('./utils'));
addpath(genpath('./forward'));

load data/synthetic_models.mat
load data/d_obs.mat


nd = length(d_obs);
nlevel = 0.1;
Cd_diag = (nlevel*d_obs).^2;
Cd = diag(Cd_diag);
inflateStd = sqrt(abs(Cd_diag));

nx = 400; nz = 40; ne = 500; N = 3;
rm = [40, 5];
nm = prod(2*rm);

niters = 4;
alphas = [9.33, 7.0, 4.0, 2.0];

pred_prod = zeros(nd, ne, niters+1);

for iter = 1:niters
    alpha = alphas(iter);
    alphaSqrt = sqrt(alpha);
    
    d_pred_ens = zeros(nd, ne);
    
    parfor iii = 1:ne
        fprintf('iteration: %d, realization: %d\n', iter, iii);
        d_pred = forward_func(EC_prior(:, :, iii), MS_prior(:, :, iii));
        d_pred_ens(:, iii) = d_pred;
    end
    
    pred_prod(:, :, iter) = d_pred_ens;    
    
    [B_EC, Qlist_EC] = rtencomp(double(EC_prior), [100, 10], 10, 2);
    [C_EC, U_EC] = hooi(B_EC, rm, 'hosvd', 4);
    [B_MS, Qlist_MS] = rtencomp(double(MS_prior), [100, 10], 10, 2);
    [C_MS, U_MS] = hooi(B_MS, rm, 'hosvd', 4);
    
    priors_rtd = [reshape(C_EC, [], ne); reshape(C_MS, [], ne)];
    
    model_mu = mean(priors_rtd, 2);
    d_pred_ens_mu = mean(d_pred_ens, 2);
    Cmd = 1/(ne-1) * (priors_rtd - model_mu) * (d_pred_ens - d_pred_ens_mu)';
    Cdd = 1/(ne-1) * (d_pred_ens - d_pred_ens_mu) * (d_pred_ens - d_pred_ens_mu)';
    
    duc = d_true + alphaSqrt*inflateStd .* randn(nd, ne);
    K = Cmd * pinv(Cdd + alpha*Cd);
    post_rtd = priors_rtd + K*(duc - d_pred_ens);
    C_post_EC = post_rtd(1:end/2, :);
    C_post_MS = post_rtd(1+end/2:end, :);
    
    C_post_MS = reshape(C_post_MS, rm(1), rm(2), ne);
    C_post_EC = reshape(C_post_EC, rm(1), rm(2), ne);
    
    B_post_MS = C_post_MS;
    B_post_EC = C_post_EC;
    
    for i = 1:N-1
        B_post_EC = ttm(B_post_EC, U_EC{i}, i);
        B_post_MS = ttm(B_post_MS, U_MS{i}, i);
    end
    
    post_EC = B_post_EC;
    post_MS = B_post_MS;
    
    for i = 1:N-1
        post_EC = ttm(post_EC, Qlist_EC{i}, i);
        post_MS = ttm(post_MS, Qlist_MS{i}, i);
    end
    
    EC_prior = post_EC;
    MS_prior = post_MS;
    
end

parfor iii = 1:ne
    fprintf('iteration: %d, realization: %d\n', iter+1, iii);
    d_pred = forward_func(EC_prior(:, :, iii), MS_prior(:, :, iii));
    d_pred_ens(:, iii) = d_pred;
end

pred_prod(:, :, iter+1) = d_pred_ens;

save data/results post_EC post_MS pred_prod