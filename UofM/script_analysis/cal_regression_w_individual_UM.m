% clear all;
% user_addpath(false,false);

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
load([save_path '/UM_info_ch128.mat']);

% idx = 4;
% filter = 'alpha';
% filter = 4;
band_name = ['band_[' num2str(bands(filter,1)) '-' ...
    num2str(bands(filter,2)) ']'];
% state = 'EO';
% state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'BS', 'DS', 'DA', 'ROC'};
% st = find(matches(state_list, state));

load('references/comb_centroids_20240311/combined_centroids_20240311.mat');
mask = C_comb.EOEC(:,topo_idx_comb);

% load(['../UM_data/data/20240402_wo_reref/topo_vector_st_' ...
%     state '_' subjname{idx} '.mat']);

save_file_path = ['../results_yjp/movie_rel_phase/' band_name '_20240723/100ms/'];
% save_file_path = ['../results_yjp/movie_rel_phase/' band_name '_20240723/20ms/'];
if ~exist([save_file_path 'regression/'], "dir")
    mkdir([save_file_path 'regression/']);
end
load([save_file_path 'topo_vector/topo_vector_st_' state '_' subjname{idx} '.mat']);

[beta, epsilon, readme] = cal_regression(mask, topo_vector);

save([save_file_path 'regression/regr_st_' state '_' subjname{idx} '.mat'], ...
    'beta', 'epsilon', 'readme', '-v7.3');

% ---------- with mae of epsilon -------------------------
if ~exist([save_file_path 'regression_mae/'], "dir")
    mkdir([save_file_path 'regression_mae/']);
end
load([save_file_path 'topo_vector/topo_vector_st_' state '_' subjname{idx} '.mat']);

[beta, epsilon, readme] = cal_regression(mask, topo_vector);
epsilon = mean(abs(epsilon),1);

save([save_file_path 'regression_mae/regr_st_' state '_' subjname{idx} '.mat'], ...
    'beta', 'epsilon', 'readme', '-v7.3');

