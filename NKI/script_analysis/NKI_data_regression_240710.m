% clear, close all;
% user_addpath(false,false);
% 
% sub = 1;
task_name = 'checkeron';
d_type = 'raw';

load('../../mfiles/references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', 'topo_idx_comb', 'topo_idx_UM2comb');
mask = C_comb.EOEC;
mask = mask(:,topo_idx_comb);

% checkeroff and checkerout cases
% without re-referencing
% load(sprintf('../preproc_%s/sub-%02d_ses-01_task-%s_eeg_topo_vector.mat', ...
%     task_name,sub,task_name));

% checkeron case
% with averaging re-referencing
% topo_vector = load(['../data/checkerboard/topo_vector_CB_subj_' num2str(sub) '.mat']).topo_vector1;
% topo_vector = topo_vector(:,topo_idx_UM2comb);

% checkeron case
% without re-referencing
load(sprintf('../data/topo_20ms/topo_100ms_wo_reref_cb/sub-%02d ses-01 task-checker eeg GAPArm 100ms_topo_wo_reref.mat', ...
    sub), 'topo');
nT = length(topo);
topo_vector = nan([nT, length(topo_idx_comb)]);
for t = 1:nT
    topo_vector(t,:) = topo{t}(topo_idx_comb);
end

K = 4;
unit = 0.1; % 100ms
n_nulls = 1000;

idx_EC = [(1:200) (1:200)+400 (1:200)+800 ...
    (1:200)+1200 (1:200)+1600];
idx_CB = [(1:200)+200 (1:200)+600 (1:200)+1000 ...
    (1:200)+1400 (1:200)+1800];

if ~exist(sprintf('../regression_%s/',task_name), "dir")
    mkdir(sprintf('../regression_%s/',task_name));
end

if ~exist(sprintf('../regression_%s/EC',task_name), "dir")
    mkdir(sprintf('../regression_%s/EC',task_name));
end
[beta, epsilon, readme] = cal_regression(mask, topo_vector(idx_EC,:));
[IDX, w_prop] = cal_regression_clustering(beta, K, unit, n_nulls);
prop = cal_transition_prop_v2(IDX, K, unit);
save(sprintf('../regression_%s/EC/sub-%02d_ses-01_task-%s_eeg_regression_st_EC.mat', ...
    task_name,sub,task_name), ...
    'IDX', 'beta', 'epsilon', 'readme', ...
    'prop', 'w_prop', 'n_nulls', 'unit', 'K', '-v7.3');

if ~exist(sprintf('../regression_%s/CB',task_name), "dir")
    mkdir(sprintf('../regression_%s/CB',task_name));
end
[beta, epsilon, readme] = cal_regression(mask, topo_vector(idx_CB,:));
[IDX, w_prop] = cal_regression_clustering(beta, K, unit, n_nulls);
prop = cal_transition_prop_v2(IDX, K, unit);
save(sprintf('../regression_%s/CB/sub-%02d_ses-01_task-%s_eeg_regression_st_CB.mat', ...
    task_name,sub,task_name), ...
    'IDX', 'beta', 'epsilon', 'readme', ...
    'prop', 'w_prop', 'n_nulls', 'unit', 'K', '-v7.3');