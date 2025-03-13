% clear;

load('../results_yjp/UM_info_AE.mat'); 
% state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'DA', 'ROC'};
load('references/topo_vector_idx_v4.mat');
topo_vector_idx = topo_vector_idx_v4;
topo_whole = [];
% idx = 3;
% filter = 'alpha';
% state = state_list{1};
% state = 'DS';
if matches(state,'LOC')
    length_t = 9000; % 1800;
else
    length_t = 15000; % 3000;
end

% tmp = topo{1};
% topo_vector_idx_v3 = find(~isnan(tmp));
% [topo_idx_x,topo_idx_y] = find(~isnan(tmp));
topo_vector = zeros(length_t,length(topo_vector_idx));
for idx = 1:n_sub_ae
% for idx = [3 5 6 7 8 9] % Burst suppression
    load(['../results_yjp/movie_rel_phase/band_alpha_whole/' subjname_ae{idx} ...
        '/' subjname_ae{idx} '_st_' state '_band_alpha_tm10_tw10_sm20_relp.mat']);
    for i = 1:length_t
        tmp = topo{i};
        topo_vector(i,:) = tmp(topo_vector_idx);
    end
    topo_whole = cat(1,topo_whole,topo_vector);
end

%%
K=4;
test_num = 10;

% k-mean clustering core
tic
disp( [ 'k=' num2str(i) ] )
[IDX,C,SUMD,D]=kmeans(topo_whole,K, 'distance', 'sqeuclidean','Replicates',test_num,'Display','final','emptyaction','drop');   
toc

%%

mode_series = zeros(n_sub_ae,length_t);
for idx = 1:n_sub_ae
    mode_series(idx,:) = IDX((idx-1)*length_t+1:idx*length_t);
end

mode_vector = C;
mode_topo   = centroid_K;
mode_dist   = D;
subj_name   = subjname_ae;

save(['../results_yjp/movie_rel_phase/mode_dynamics_st_' state '_K_' num2str(K) '.mat'], ...
    'mode_series', 'mode_vector', 'mode_topo', 'mode_dist', ...
    'subj_name', 'chan_coord_xy', 'shiftpreset', 'smooth', ...
    'time_window', 'time_window', 'time_all', 'time_pt');