% clear all;
user_addpath(false,false);

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
load([save_path '/UM_info_ch128.mat']);

% idx = 4;
% filter = 'alpha';
% filter = 4;
band_name = ['band_[' num2str(bands(filter,1)) '-' ...
    num2str(bands(filter,2)) ']'];
% unit = 0.1; % 100ms
% unit = 0.02; % 20ms
if ismember(filter, [2, 3]) % delta(1-4), theta(4-8)
    unit = 0.1;
elseif ismember(filter, [4, 5, 6, 7]) % alpha, low/high beta, gamma
    unit = 0.02;
end

% state = 'EO';
% state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'BS', 'DS', 'DA', 'ROC'};
% st = find(matches(state_list, state));

% load(['../UM_data/regression/20240402_wo_reref/regr_st_' ...
%     state '_' subjname{idx} '.mat']);
save_file_path = ['../results_yjp/movie_rel_phase/' band_name '_20240723/' num2str(unit*1e3) 'ms/'];
if ~exist([save_file_path 'clustering_w_1min/'], "dir")
    mkdir([save_file_path 'clustering_w_1min/']);
end
load([save_file_path 'regression/regr_st_' state '_' subjname{idx} '.mat']);

K = 4;
n_nulls = 1000;

slice_len = fix(60/unit); % 1min:60sec
slice_num = size(beta,2)/slice_len;

IDX = struct();
prop = struct();
w_prop = struct();

for sl = 1:slice_num
    [IDX.(sprintf('slice_%d',sl)), w_prop.(sprintf('slice_%d',sl))] = cal_regression_clustering( ...
        beta(:,1+(sl-1)*slice_len:sl*slice_len), K, unit, n_nulls);
    prop.(sprintf('slice_%d',sl)) = cal_transition_prop_v2(IDX.(sprintf('slice_%d',sl)), K, unit);
end

save([save_file_path 'clustering_w_1min/regr_st_' state '_' subjname{idx} '_n_nulls_' num2str(n_nulls) '.mat'], ...
    'IDX', 'w_prop', 'prop', ...
    'slice_num', 'slice_len', 'unit', 'K', 'n_nulls', '-v7.3');

