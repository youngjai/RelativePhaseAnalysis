% clear, close all;
% user_addpath(false, false);
% 
% sub = 1;
% task_type = 'checkeron';

load('../../UM_data/pca/topo_20240723/pca_st_EC_all.mat', 'coeff');
load('../../mfiles/references/comb_centroids_20240311/combined_centroids_20240311.mat');
mask = C_comb.EOEC(:,topo_idx_comb);

if ~matches(task_type(end-1:end), 'on')
    % checkeroff and checkerout
    load(sprintf('../preproc_%s/sub-%02d_ses-01_task-%s_eeg_topo_vector.mat', ...
        task_type, sub, task_type), 'topo_vector');
else
    % checkeron
    load(sprintf('../data/topo_20ms/topo_100ms_wo_reref_cb/sub-%02d ses-01 task-checker eeg GAPArm 100ms_topo_wo_reref.mat', ...
        sub), 'topo');
    topo_vector = nan([length(topo),length(topo_idx_comb)]);
    for t = 1:length(topo)
        topo_vector(t,:) = topo{t}(topo_idx_comb);
    end
end

idx_EC = [(1:200) (1:200)+400 (1:200)+800 ...
    (1:200)+1200 (1:200)+1600];
idx_CB = [(1:200)+200 (1:200)+600 (1:200)+1000 ...
    (1:200)+1400 (1:200)+1800];

tic;
save_path = ['../PCA_' task_type '/'];
if ~exist(save_path, "dir")
    mkdir(save_path);
end
state = 'CB';
if ~exist([save_path state '/'], "dir")
    mkdir([save_path state '/']);
end
score = topo_vector(idx_CB,:)*coeff;
latent = var(score,0,1)';
explained = 100*latent/sum(latent);
save([save_path state ...
    sprintf('/sub-%02d_ses-01_task-%s_eeg_pca_st_%s.mat', ...
    sub, task_type, state)], ...
    'score', 'latent', 'explained', 'idx_CB', '-v7.3');

state = 'EC';
if ~exist([save_path state '/'], "dir")
    mkdir([save_path state '/']);
end
score = topo_vector(idx_EC,:)*coeff;
latent = var(score,0,1)';
explained = 100*latent/sum(latent);
save([save_path state ...
    sprintf('/sub-%02d_ses-01_task-%s_eeg_pca_st_%s.mat', ...
    sub, task_type, state)], ...
    'score', 'latent', 'explained', 'idx_EC', '-v7.3');

toc;
