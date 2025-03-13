%% pca mask and indexing (send to Younghwa, for SI)
% clear, close all;
% user_addpath(false, false); % set path function (you don't need it.)

% load topo_vector and concatenate it
load('../results_yjp/UM_info_ch128.mat', 'subjname', 'n_sub');
load_path = ['../results_yjp/movie_rel_phase/' band_name '_20240723'];
% resol = 20;
unit = resol/1e3;

topo_vector_all = [];
for idx = 1:n_sub
    try
        disp(idx)
        tic;
        topo_vector_all = cat(1,topo_vector_all, ...
            load(sprintf('%s/%dms/topo_vector/topo_vector_st_EC_%s.mat', ...
            load_path,resol,subjname{idx})).topo_vector);
        toc;
    end
end

% %%
% perform pca without centerization
tic;
% shape of topo_vector must be 'observed points' by 'variables'
[coeff, score, latent, tsquared, explained, mu] = pca(topo_vector_all,'Centered',false);
toc;

% save pca results
save(sprintf('%s/%dms/pca_result_of_EC.mat',load_path,resol), ...
    'coeff','score','latent','tsquared','explained','mu','-v7.3');
save(sprintf('%s/%dms/pca_result_of_EC_explained.mat',load_path,resol), ...
    'explained','-v7.3');

% %%
% clear, close all;
% user_addpath(false, false); % set path function (you don't need it.)

% load(sprintf('%s/%dms/pca_result_of_EC.mat',load_path,resol));
load('references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'topo_idx_comb', 'C_comb');
load('references/topo_coords_ch96.mat');

n_pc = 8;

eig_vec = zeros([n_pc,size(score,2)]); 
for pc = 1:n_pc
    eig_vec(pc,pc) = 1; 
end
eig_map = eig_vec*coeff';

% % make mask to conduct regression using 'cal_regression.m' function
% load('references/comb_centroids_20240311/combined_centroids_20240311.mat','C_comb');
check_corr = corr(C_comb.EOEC(1:2,topo_idx_comb)',eig_map(1:2,:)');
% step 1: change index (pc1:front-to-back, pc2:left-to-right)
if check_corr(1)<check_corr(2)
    eig_map(1:2,:) = eig_map([2 1],:);
    check_corr = corr(C_comb.EOEC(1:2,topo_idx_comb)',eig_map(1:2,:)');
    explained(1:2) = explained([2 1]);
end
% step 2: reverse the map (e.g. back-to-front -> front-to-back)
eig_map(1:2,:) = eig_map(1:2,:).*sign(diag(check_corr));

% %%
% draw and save principal components
% user_addpath(false, false); % set path function (you don't need it.)
% load('references/comb_centroids_20240311/combined_centroids_20240311.mat', 'topo_idx_comb');
% load('references/topo_coords_ch96.mat');

fig = figure(1);
fig.Position = [100 100 800 500];
clf;
for pc = 1:n_pc
    subplot_tight(2,n_pc/2,pc);
    topo = nan(200);
    topo(topo_idx_comb) = eig_map(pc,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord);
    caxis([-max(abs(caxis)) max(abs(caxis))]);
    % if caxis doesn't work, try bellow command:
    % clim([-max(abs(clim)) max(abs(clim))]);
    title(sprintf('PC %d (%.2f%%)',pc,explained(pc)), ...
        'FontSize',10, 'FontWeight','normal');
    colormap('jet');
    colorbar;
end
sgtitle(band_name, 'Interpreter','none');

% save pc figure
exportgraphics(fig, sprintf('%s/%dms/pca_result_of_EC.png',load_path,resol));

% save pc mask
mask = [eig_map(1:2,:); -eig_map([2 1],:)];
save(sprintf('%s/%dms/pca_mask_result_of_EC.mat',load_path,resol), ...
    'mask', '-v7.3');

% % %% to get result of each subject...
% % clear, close all;
% 
% % load(sprintf('%s/%dms/pca_mask_result_of_EC.mat',load_path,resol));
% 
% [beta, epsilon, readme] = cal_regression(mask, topo_vector);
% epsilon_mae = mean(abs(epsilon),1,'omitnan');
% 
% K = 4;
% % unit = 0.1; % 100ms
% % unit = 0.02; % 20ms
% n_nulls = 1000;
% % calculate weighted statistics
% [IDX, w_prop] = cal_regression_clustering(beta, K, unit, n_nulls);
% % calculate unweighted statistics
% prop = cal_transition_prop_v2(IDX, K, unit);
% 
% save(sprintf('%s/%dms/pca_mask_result_of_EC.mat',load_path,resol), ...
%     'IDX', 'w_prop', 'prop', '-v7.3');