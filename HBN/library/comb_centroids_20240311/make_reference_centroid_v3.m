% to make reference centroids with EO and EC
% EO : 70(HBN, 16.8*5=84s, 5544s), 18(UofM, 300s, 5100s) = 88
% EC : 70(HBN, 33.6*5=168s, 11088s), 18(UofM, 300s, 5100s) = 88
% EOEC : 88 + 88 = 176

% 12, 54줄에서 n_HBN 숫자 바꾸고,
% 24, 26, 34, 36, 49, 56, 57, 76, 115줄 경로 바꾸세요~!

%% combined centroids (HBN, UofM) EO and EC
clear all; close all;

n_UM = 18;
n_HBN = 66;
n_total = n_UM + n_HBN;

K = 4;

state_list = {'EC','EO'};
n_state = length(state_list);

for st = 1:n_state
    state = state_list{st};
    
    topo_idx_UM = load('references/topo_vector_idx_v4.mat', ...
        'topo_vector_idx_v4').topo_vector_idx_v4;
    C_UM = load(['../UM_data/k_means/20240223_wo_reref/k_means_' ...
        state '_all.mat'], 'C').C;
    C_UM_topo = nan([4,200,200]);
    
    for i = 1:length(C_UM)
        C_UM_topo(:,topo_idx_UM(i)) = C_UM(:,i);
    end
    
    topo_idx_HBN = load('/youngjai/Dropbox/Moonbrainlab/HBN/centroid_topo_ref/topo_sort_ref.mat', ...
        'topo_sort_ref').topo_sort_ref;
    C_HBN = load(['/youngjai/Dropbox/Moonbrainlab/HBN/centroid_topo_ref/HBN_control_' ...
        lower(state) '.mat'], 'C').C;
    C_HBN_topo = nan([4,200,200]);
    
    for i = 1:length(C_HBN)
        C_HBN_topo(:,topo_idx_HBN(i)) = C_HBN(:,i);
    end

    C_comb.(state) = (n_UM/n_total).*C_UM_topo + (n_HBN/n_total).*C_HBN_topo;
end

C_comb.EOEC = 0.5.*(C_comb.EO + C_comb.EC);

[topo_idx_comb, topo_idx_UM2comb, topo_idx_HBN2comb] = ...
    intersect(topo_idx_UM, topo_idx_HBN);

save('references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', 'topo_idx_UM', 'topo_idx_HBN', ...
    'topo_idx_comb', 'topo_idx_UM2comb', 'topo_idx_HBN2comb', ...
    'n_UM', 'n_HBN', 'n_total', '-v7.3');

%% draw combined centroids
clear all; close all;

load('references/comb_centroids_20240311/combined_centroids_20240311.mat');
load('references/topo_coords_ch96.mat');

K = 4;

state_list = {'EO', 'EC', 'EOEC'};

fig = figure(1);
fig.Position = [100 100 1500 450];
for st = 1:length(state_list)
    state = state_list{st};
    clf;
    sgtitle(['Combined centroids of ' state ', clim([-0.6 0.6])'], ...
        'FontSize',20, 'FontWeight','bold');
    for i = 1:K
        subplot(1,4,i);
        topoplot_figure(squeeze(C_comb.(state)(i,:,:)), ...
            borderCoords, xx, yy, Coord, 'scatter',1);
        clim([-0.6 0.6]);
    end
    exportgraphics(fig, ['references/comb_centroids_20240311/combined_centroids_20240311_' state '.png']);
end

% correlation between them

corr_mat = zeros([K,K,3]);
pval_mat = zeros([K,K,3]);

state_list1 = cell([3,2]);
state_list1(:,1) = {'EO', 'EO', 'EC'};
state_list1(:,2) = {'EC', 'EOEC', 'EOEC'};

for st = 1:3
    st1 = state_list1{st,1};
    st2 = state_list1{st,2};
    for k1 = 1:K
        for k2 = 1:K
            x = C_comb.(st1)(k1,isfinite(C_comb.(st1)(k1,:)))';
            y = C_comb.(st2)(k2,isfinite(C_comb.(st2)(k2,:)))';
            [corr_mat(k1,k2,st), pval_mat(k1,k2,st)] = corr(x,y);
        end
    end
end

fig = figure(1);
fig.Position = [100 100 1500 400];
clf;
for st = 1:length(state_list1)
    st1 = state_list1{st,1};
    st2 = state_list1{st,2};
    subplot(1,length(state_list1),st);
    h = heatmap(corr_mat(:,:,st));
    h.XLabel = st2;
    h.YLabel = st1;
    h.ColorLimits = [-1, 1];
    h.Colormap = turbo;
    h.CellLabelFormat = '%7.4f';
    set(gca, 'FontSize',15);
end
exportgraphics(fig, ['references/comb_centroids_20240311/combined_centroids_20240311_corr_mat.png']);

