%% load file and setup

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
% load([save_path '/UM_info_AE.mat']);
% subjname = subjname_ae;
load([save_path '/UM_info.mat']);

load('references/topo_vector_idx_v4.mat');
load('references/topo_coords_ch96.mat');

% Define file name and load the file.
% idx = 1;
% state = 'EO';
% filter = 'alpha';
time_moving = 50;
time_window = 50;

smooth = 20;
% vectorize data
load([save_path '/movie_rel_phase/band_alpha_re_whole_w_EGIchanlocs/' subjname{idx} ...
    '/' subjname{idx} '_st_' state '_band_' filter '_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo.mat']);

topo_size = length(topo);
topo_vector = zeros(topo_size,length(topo_idx_x));

% topo_vector: time point x frame vector
for i=1:topo_size
    T = topo{i};
    topo_vector(i,:) = T(topo_vector_idx_v4);
end

shiftpreset = 0;
% 0 if it is human
% 1 if it is monkey

%% K-mean clustering
K=4;
% sort 전 좌 후 뒤
% mask = load(['references/mask_st_EO_K_4_band_alpha_tm' num2str(time_moving) ...
%     '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo.mat']);
mask = load('../UM_data/k_means/k_means_EC_all.mat'] 'C').C;
[IDX,C,SUMD,D] = kmeans_with_mask_v3(mask, topo_vector, 'euclidean', 1);
% [IDX,C,SUMD,D] = kmeans_with_mask_v2(mask.C, topo_vector, 'pearson');

file_path = [save_path '/k_means'];
if ~exist(file_path, "dir")
    mkdir(file_path);
end
file_path = [save_path '/k_means/band_' filter '_re_whole_w_EGIchanlocs_w_EC_w_noto'];
if ~exist(file_path, "dir")
    mkdir(file_path);
end
if ~exist([file_path '/' subjname{idx}], "dir")
    mkdir([file_path '/' subjname{idx}]);
end

save([file_path '/' subjname{idx} ...
    '/' subjname{idx} '_st_' state '_band_' filter '_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans.mat'], ...
    'IDX', 'C', 'SUMD', 'D', 'time_moving', 'time_window', 'smooth','-v7.3');

%%
if ~exist([file_path '/figrues'], "dir")
    mkdir([file_path '/figures']);
end
if ~exist([file_path '/figures/' subjname{idx}], "dir")
    mkdir([file_path '/figures/' subjname{idx}]);
end


fig = figure('Visible','off');
fig.Position(3:4) = [600,500];
clf;
histogram(IDX,'normalization','probability')
xlabel('$K_i$', 'FontSize',20, 'Interpreter','latex');
ylabel('${\rm PDF}$', 'FontSize',20, 'Interpreter','latex');
exportgraphics(fig, [file_path '/figures/' subjname{idx} ...
    '/' subjname{idx} '_st_' state '_band_' filter '_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans_modes.png']);


fig.Position(3:4) = [1200,400];
clf;
plot(time_all, IDX,'.-')
xlabel('$t ({\rm second})$', 'FontSize',20, 'Interpreter','latex');
ylabel('$K_i$', 'FontSize',20, 'Interpreter','latex');
exportgraphics(fig, [file_path '/figures/' subjname{idx} ...
    '/' subjname{idx} '_st_' state '_band_' filter '_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans_series.png']);
