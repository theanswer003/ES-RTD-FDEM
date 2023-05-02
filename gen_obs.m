clear; clc; close all;
addpath(genpath('./forward'));
load ./data/ObservedData_1Model.mat

d_true = [FWD_IP_RC'; FWD_QP_RC'];
d_true = d_true';
d_true = d_true(:);
nlevel = 0.05;
d_obs = d_true + nlevel * d_true .* randn(size(d_true));

save data/d_obs d_obs d_true