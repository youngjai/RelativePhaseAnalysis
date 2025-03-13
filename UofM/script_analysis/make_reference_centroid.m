% to make reference centroids with EO and EC
% EO : 66(HBN, 16.8*5=84s, 5544s), 17(UofM, 300s, 5100s), 
% 17(NKI, 2.1*288=604.8s, 10281.6s) = 100
% EC : 66(HBN, 33.6*5=168s, 11088s), 18(UofM, 300s, 5400s) = 84
% EOEC : 100 + 84 = 184

%% combined centroids (HBN, UofM, NKI) EO
clear all; close all;

load('references/topo_vector_idx_v4.mat');
load('../HBN_data/data/topo_vector_idx_HBN.mat', 'topo_vector_idx_HBN');
[~,reidx_HBN] = intersect(topo_vector_idx_HBN,topo_vector_idx_v4);

K = 4;
comb_C = zeros([K,length(topo_vector_idx_v4)]);
n_total = 0;

% HBN
n_subj = 66;
n_subj1 = 0;
for sub = 1:n_subj
    load(['../HBN_data/k_means/k_means_EO_subj_' num2str(sub) '.mat'], 'C');
    C = C(:,reidx_HBN);
    comb_C = comb_C + C;
    n_subj1 = n_subj1 + 1;
end
disp(['HBN completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% UofM
load('../results_yjp/UM_info.mat', 'subjname', 'labels');
n_subj = length(subjname);
n_subj1 = 0;
for sub = 1:n_subj
    try
        load(['../results_yjp/k_means/band_alpha_re_whole/' subjname{sub} ...
            '/' subjname{sub} '_st_EO_band_alpha_tm50_tw50_sm20_topo_kmeans.mat'], 'C');
        disp([subjname{sub} ' ' labels{sub}]);
        comb_C = comb_C + C;
        n_subj1 = n_subj1 + 1;
    end
end
disp(['UofM completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% NKI
reidx_NKI = [1 2 4 6 7 8 9 12 13 14 15 16 17 18 20 21 22];
n_subj = length(reidx_NKI);
n_subj1 = 0;
for sub = 1:n_subj
    load(['../NKI_data/k_means/k_means_subj_' num2str(reidx_NKI(sub)) '.mat'], 'C');
    comb_C = comb_C + C;
    n_subj1 = n_subj1 + 1;
end
disp(['NKI completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;
comb_C = comb_C./n_total;
disp(['All completed, N:' num2str(n_total)]);

save('references/combined_centroids_EO.mat', 'comb_C', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');

%% combined centroids (HBN, UofM, NKI) EC
clear all; close all;

load('references/topo_vector_idx_v4.mat');
load('../HBN_data/data/topo_vector_idx_HBN.mat', 'topo_vector_idx_HBN');
[~,reidx_HBN] = intersect(topo_vector_idx_HBN,topo_vector_idx_v4);

K = 4;
comb_C = zeros([K,length(topo_vector_idx_v4)]);
n_total = 0;

% HBN
n_subj = 66;
n_subj1 = 0;
for sub = 1:n_subj
    load(['../HBN_data/k_means/k_means_EC_subj_' num2str(sub) '.mat'], 'C');
    C = C(:,reidx_HBN);
    comb_C = comb_C + C;
    n_subj1 = n_subj1 + 1;
end
disp(['HBN completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% UofM
load('../results_yjp/UM_info.mat', 'subjname', 'labels');
n_subj = length(subjname);
n_subj1 = 0;
for sub = 1:n_subj
    try
        load(['../results_yjp/k_means/band_alpha_re_whole/' subjname{sub} ...
            '/' subjname{sub} '_st_EC_band_alpha_tm50_tw50_sm20_topo_kmeans.mat'], 'C');
        disp([subjname{sub} ' ' labels{sub}]);
        comb_C = comb_C + C;
        n_subj1 = n_subj1 + 1;
    end
end
disp(['UofM completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;
comb_C = comb_C./n_total;
disp(['All completed, N:' num2str(n_total)]);

save('references/combined_centroids_EC.mat', 'comb_C', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');


%% combined centroids (HBN, UofM, NKI) EOEC
clear all; close all;

load('references/combined_centroids_EO.mat');
EC = load('references/combined_centroids_EC.mat', 'comb_C', 'n_total');

comb_C = 0.5*(comb_C+EC.comb_C);
n_total = n_total + EC.n_total;
save('references/combined_centroids_EOEC.mat', 'comb_C', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');

%% draw combined centroids
clear all; close all;

% state = 'EO';
state_list = {'EO', 'EC', 'EOEC'};

load('references/topo_vector_idx_v4.mat');
load('references/topo_coords_ch95.mat');

for st = 1:length(state_list)
    state = state_list{st};
    load(['references/combined_centroids_' state '.mat']);
    K = size(comb_C,1);
    centroid_K = nan(size(xx,1),size(xx,2),K);
    
    for j=1:K
        for i=1:length(topo_idx_x)
            centroid_K(topo_idx_x(i),topo_idx_y(i),j) = comb_C(j,i); 
        end
    end
    
    fig = figure(1);
    fig.Position = [100 100 1400 450];
    clf;
    sgtitle(['Combined centroids with ' state ', color range [-0.5, 0.5]'], 'FontSize',20);
    for k=1:K
        subplot(1,K,k);
        topoplot_figure(centroid_K(:,:,k), borderCoords, xx, yy, Coord, 'scatter', 1);
        title(['Mode ' num2str(k)], 'FontSize',15);
        clim([-0.5 0.5]);
    end
    exportgraphics(fig, ['references/combined_centroids_' state '.png']);
end

%% correlation between centroids 1
clear all; close all;

state_list = {'EO', 'EC', 'EOEC'};
C = struct;

for st = 1:length(state_list)
    state = state_list{st};
    C.(state) = load(['references/combined_centroids_' state '.mat'], 'comb_C').comb_C;
end

%% correlation between centroids 2
state1 = 'EO'; state2 = 'EC';
% state1 = 'EC'; state2 = 'EOEC';
% state1 = 'EOEC'; state2 = 'EO';
fig = figure(1);
fig.Position = [100 100 600 600];
corr_mat = corr(C.(state1)', C.(state2)');
h = heatmap(corr_mat);
h.XLabel = state2;
h.YLabel = state1;
h.FontSize = 20;
h.ColorLimits = [-1, 1];
h.CellLabelFormat = '%3.4f';
h.Colormap = cat(1,bone,flipud(bone));
h.Title = sprintf('mean corr. coef.: %3.4f', mean(diag(corr_mat)));
exportgraphics(fig, ['references/combined_centroids_correlation_' ...
    state1 '_and_' state2 '.png']);


%% combined principal components (HBN, UofM, NKI) EO
clear all; close all;

load('references/topo_vector_idx_v4.mat');
load('../HBN_data/data/topo_vector_idx_HBN.mat', 'topo_vector_idx_HBN');
[~,reidx_HBN] = intersect(topo_vector_idx_HBN,topo_vector_idx_v4);

K = 4;
n_pnts = 50;
comb_coeff = zeros([length(topo_vector_idx_v4),n_pnts]);
n_total = 0;

% HBN
n_subj = 66;
n_subj1 = 0;
for sub = 1:n_subj
    load(['../HBN_data/pca/pca_EO_subj_' num2str(sub) '.mat'], 'coeff');
    coeff = coeff(reidx_HBN,1:n_pnts);
    comb_coeff = comb_coeff + coeff;
    n_subj1 = n_subj1 + 1;
end
disp(['HBN completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% UofM
load('../results_yjp/UM_info.mat', 'subjname', 'labels');
n_subj = length(subjname);
n_subj1 = 0;
for sub = 1:n_subj
    try
        load(['../results_yjp/pca/' subjname{sub} ...
            '/' subjname{sub} '_st_EO_band_alpha_tm50_tw50_sm20_topo_kmeans_pca.mat'], 'coeff');
        disp([subjname{sub} ' ' labels{sub}]);
        coeff = coeff(:,1:n_pnts);
        comb_coeff = comb_coeff + coeff;
        n_subj1 = n_subj1 + 1;
    end
end
disp(['UofM completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% NKI
reidx_NKI = [1 2 4 6 7 8 9 12 13 14 15 16 17 18 20 21 22];
n_subj = length(reidx_NKI);
n_subj1 = 0;
for sub = 1:n_subj
    load(['../NKI_data/pca/pca_rest_subj_' num2str(reidx_NKI(sub)) '.mat'], 'coeff');
    coeff = coeff(:,1:n_pnts);
    comb_coeff = comb_coeff + coeff;
    n_subj1 = n_subj1 + 1;
end
disp(['NKI completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;
comb_coeff = comb_coeff./n_total;
disp(['All completed, N:' num2str(n_total)]);

save('references/combined_PCs_EO.mat', 'comb_coeff', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');

%% combined principal components (HBN, UofM, NKI) EC
clear all; close all;

load('references/topo_vector_idx_v4.mat');
load('../HBN_data/data/topo_vector_idx_HBN.mat', 'topo_vector_idx_HBN');
[~,reidx_HBN] = intersect(topo_vector_idx_HBN,topo_vector_idx_v4);

K = 4;
n_pnts = 50;
comb_coeff = zeros([length(topo_vector_idx_v4),n_pnts]);
n_total = 0;

% HBN
n_subj = 66;
n_subj1 = 0;
for sub = 1:n_subj
    load(['../HBN_data/pca/pca_EC_subj_' num2str(sub) '.mat'], 'coeff');
    coeff = coeff(reidx_HBN,1:n_pnts);
    comb_coeff = comb_coeff + coeff;
    n_subj1 = n_subj1 + 1;
end
disp(['HBN completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;

% UofM
load('../results_yjp/UM_info.mat', 'subjname', 'labels');
n_subj = length(subjname);
n_subj1 = 0;
for sub = 1:n_subj
    try
        load(['../results_yjp/pca/' subjname{sub} ...
            '/' subjname{sub} '_st_EC_band_alpha_tm50_tw50_sm20_topo_kmeans_pca.mat'], 'coeff');
        disp([subjname{sub} ' ' labels{sub}]);
        coeff = coeff(:,1:n_pnts);
        comb_coeff = comb_coeff + coeff;
        n_subj1 = n_subj1 + 1;
    end
end
disp(['UofM completed, N:' num2str(n_subj1)]);
n_total = n_total + n_subj1;
comb_coeff = comb_coeff./n_total;
disp(['All completed, N:' num2str(n_total)]);

save('references/combined_PCs_EC.mat', 'comb_coeff', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');

%% combined principal components (HBN, UofM, NKI) EOEC
clear all; close all;

load('references/combined_PCs_EO.mat');
EC = load('references/combined_PCs_EC.mat', 'comb_coeff', 'n_total');

comb_coeff = 0.5*(comb_coeff+EC.comb_coeff);
n_total = n_total + EC.n_total;
save('references/combined_PCs_EOEC.mat', 'comb_coeff', 'n_total', ...
    'topo_vector_idx_v4', 'topo_vector_idx_HBN', 'reidx_HBN');

%% draw combined PC figures
L=4;
load('references/topo_vector_idx_v4.mat');
load('references/topo_coords_ch95.mat');

state_list = {'EO', 'EC', 'EOEC'};
% state = 'EC';
for st = 1:length(state_list)
    state = state_list{st};

    load(['references/combined_PCs_' state '.mat']);
    
    pca_L = nan(size(xx,1),size(xx,2),L);
    
    for j=1:L
        for i=1:length(topo_idx_x)
            pca_L(topo_idx_x(i),topo_idx_y(i),j) = comb_coeff(i,j); 
        end
    end
    
    fig = figure(507);
    fig.Position = [100 100 1400 450];
    clf;
    sgtitle(['Combined PCs with ' state ', color range [-0.005, 0.005]'], 'FontSize',20);
    for j=1:L
        subplot(1,L,j);
        topoplot_figure(pca_L(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
        title(['PC ' num2str(j)], 'FontSize',15);
        clim([-0.005 0.005]);
    end
    exportgraphics(fig, ['references/combined_PCs_' state '.png']);
end

%% Mode distribution with all dataset
clear all; close all;

K = 4;
mode_IDX = NaN([K,12,66]); % mode, state, subject
% state
% EO, EC, AI1, LOC, AI2, BS, DS, DA, ROC
state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'BS', 'DS', 'DA', 'ROC'};

% NKI data
reidx_NKI = [1 2 4 6 7 8 9 12 13 14 15 16 17 18 20 21 22];
n_subj = length(reidx_NKI);
for subj = 1:n_subj
    load(['../NKI_data/k_means/k_means_rest_subj_' num2str(reidx_NKI(subj)) ...
        '_w_comb_mask.mat'], 'IDX');
    for k = 1:K
        mode_IDX(k,1,subj) = sum(IDX==k)./length(IDX);
    end
end

% HBN data
n_subj = 66;
for st = 1:2
    state = state_list{st};
    for subj = 1:n_subj
        load(['../HBN_data/k_means/k_means_' state ...
            '_subj_' num2str(subj) '_w_comb_mask.mat'], 'IDX');
        for k = 1:K
            mode_IDX(k,st+1,subj) = sum(IDX==k)./length(IDX);
        end
    end
end

% UofM data
load('../results_yjp/UM_info.mat', 'subjname');
n_subj = length(subjname);
for st = 1:length(state_list)
    state = state_list{st};
    for subj = 1:n_subj
        try
            load(['../results_yjp/k_means/band_alpha_re_whole_w_comb/' ...
                subjname{subj} '/' subjname{subj} '_st_' state ...
                '_band_alpha_tm50_tw50_sm20_topo_kmeans.mat'], 'IDX');
            for k = 1:K
                mode_IDX(k,st+3,subj) = sum(IDX==k)./length(IDX);
            end
        end
    end
end

mode_pdf = mean(mode_IDX,3,'omitnan');
mode_se = std(mode_IDX,1,3,'omitnan')./sqrt(sum(~isnan(mode_IDX(1,:,:)),3));

state_label_list = {'EO (NKI)', ...
    'EO (HBN)', 'EC (HBN)', ...
    'EO', 'EC', 'Propofol', 'LOC', 'Isoflurane', 'Burst', 'Suppr.', 'Anesth.', 'ROC'};

save('references/all_dataset_mode_dist.mat', ...
    'mode_IDX', 'mode_pdf', 'mode_se', 'state_label_list');

%% Mode distribution with all dataset (draw)
clear all; close all;

load('references/all_dataset_mode_dist.mat');
K = 4;
c_list = [[119 175 216]; [183 143 191]; [199 128 143]; [238 208 136]; [179 206 145]];
c_list = c_list./255;

fig = figure(1);
fig.Position = [100 100 1800 500];
clf; hold on; grid on;
x = [1 2.5:3.5 5:13];
for k = 1:K
    b(k) = bar(x+14*(k-1), mode_pdf(k,:), ...
        'FaceColor',c_list(k,:), 'LineWidth',1);
    errorbar(x+14*(k-1), mode_pdf(k,:), mode_se(k,:), ...
        'k', 'LineStyle','none');
end
x1 = [1.75 4.25];
xline([x1 x1+14 x1+28 x1+42], 'LineWidth',2);
xline([14 28 42], '--', 'LineWidth',3);
xticks([x x+14 x+28 x+42]);
xticklabels([state_label_list state_label_list state_label_list state_label_list]);
ylim([0,0.5]);
set(gca, 'FontSize',15);
xlabel('Mode');
ylabel('PDF');
legend(b, {'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4'}, ...
    'NumColumns',4);

exportgraphics(fig, ['references/all_dataset_mode_dist.png']);

%% Mode distribution with NOTA with all dataset
clear all; close all;

K = 4;
mode_IDX = NaN([K+1,12,66]); % mode, state, subject
% state
% EO, EC, AI1, LOC, AI2, BS, DS, DA, ROC
state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'BS', 'DS', 'DA', 'ROC'};

% NKI data
reidx_NKI = [1 2 4 6 7 8 9 12 13 14 15 16 17 18 20 21 22];
n_subj = length(reidx_NKI);
for subj = 1:n_subj
    load(['../NKI_data/k_means/k_means_rest_subj_' num2str(reidx_NKI(subj)) ...
        '_w_comb_mask_w_NOTA.mat'], 'IDX');
    for k = 1:K+1
        mode_IDX(k,1,subj) = sum(IDX==k)./length(IDX);
    end
end

% HBN data
n_subj = 66;
for st = 1:2
    state = state_list{st};
    for subj = 1:n_subj
        load(['../HBN_data/k_means/k_means_' state ...
            '_subj_' num2str(subj) '_w_comb_mask_w_NOTA.mat'], 'IDX');
        for k = 1:K+1
            mode_IDX(k,st+1,subj) = sum(IDX==k)./length(IDX);
        end
    end
end

% UofM data
load('../results_yjp/UM_info.mat', 'subjname');
n_subj = length(subjname);
for st = 1:length(state_list)
    state = state_list{st};
    for subj = 1:n_subj
        try
            load(['../results_yjp/k_means/band_alpha_re_whole_w_comb_w_NOTA/' ...
                subjname{subj} '/' subjname{subj} '_st_' state ...
                '_band_alpha_tm50_tw50_sm20_topo_kmeans.mat'], 'IDX');
            for k = 1:K+1
                mode_IDX(k,st+3,subj) = sum(IDX==k)./length(IDX);
            end
        end
    end
end

mode_pdf = mean(mode_IDX,3,'omitnan');
mode_se = std(mode_IDX,1,3,'omitnan')./sqrt(sum(~isnan(mode_IDX(1,:,:)),3));

state_label_list = {'EO (NKI)', ...
    'EO (HBN)', 'EC (HBN)', ...
    'EO', 'EC', 'Propofol', 'LOC', 'Isoflurane', 'Burst', 'Suppr.', 'Anesth.', 'ROC'};

save('references/all_dataset_mode_dist_w_NOTA.mat', ...
    'mode_IDX', 'mode_pdf', 'mode_se', 'state_label_list');

%% Mode distribution with NOTA with all dataset (draw)
clear all; close all;

load('references/all_dataset_mode_dist_w_NOTA.mat');
K = 4;
c_list = [[119 175 216]; [183 143 191]; [199 128 143]; [238 208 136]; [179 206 145]];
c_list = c_list./255;

fig = figure(1);
fig.Position = [100 100 1800 500];
clf; hold on; grid on;
x = [1 2.5:3.5 5:13];
for k = 1:K+1
    b(k) = bar(x+14*(k-1), mode_pdf(k,:), ...
        'FaceColor',c_list(k,:), 'LineWidth',1);
    errorbar(x+14*(k-1), mode_pdf(k,:), mode_se(k,:), ...
        'k', 'LineStyle','none');
end
x1 = [1.75 4.25];
xline([x1 x1+14 x1+28 x1+42 x1+56], 'LineWidth',2);
xline([14 28 42 56], '--', 'LineWidth',3);
xticks([x x+14 x+28 x+42 x+56]);
xticklabels([state_label_list state_label_list state_label_list state_label_list state_label_list]);
ylim([0,0.5]);
set(gca, 'FontSize',15);
xlabel('Mode');
ylabel('PDF');
legend(b, {'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5'}, ...
    'NumColumns',5);

exportgraphics(fig, ['references/all_dataset_mode_dist_w_NOTA.png']);

