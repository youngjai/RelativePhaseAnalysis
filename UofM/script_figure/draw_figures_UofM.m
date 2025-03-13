%% draw combined centroid (UofM and HBN) (save)
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/' ...
    'preprocessed_wo_reref/band_[8-12]/'],dropbox_path);
load([load_path 'UM_3_wo_badchan_and_wo_reref_and_filtered']);

state = 'EC';
st = strcmp(state_list,state);
rel_p = cal_rel_phase_v3(double(data_state{st}));
chan_coord_xy = chanlocs{st};
smooth = 20;
shiftpreset = 0;
% 0 if it is human
% 1 if it is monkey
[~, xx, yy, Coord, borderCoords] = topoplot_general_test(rel_p(1,:)', [-[chan_coord_xy.Y]; [chan_coord_xy.X]]', ...
     'smooth',smooth, 'shiftpreset',shiftpreset, 'scatter',1);
save('references/topoplot_parameters.mat', ...
    'xx', 'yy', 'Coord', 'borderCoords', '-v7.3');

%% snapshot of topoplot in EC, DS in UM_8
% (A) Phase maps over time in the conscious and unconscious state
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'BS', 'DS', 'DA', 'ROC'};
load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/' ...
    'preprocessed_wo_reref/band_[8-12]/100ms/topo_vector/'],dropbox_path);
load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'topo_idx_comb');
load('references/topoplot_parameters.mat');

state1 = 'EC';
st1 = strcmp(state_list, state1);
topo_vector1 = load([load_path 'topo_vector_st_' state1 '_UM_8.mat'], ...
    'topo_vector').topo_vector;

state2 = 'DS';
st2 = strcmp(state_list, state2);
topo_vector2 = load([load_path 'topo_vector_st_' state2 '_UM_8.mat'], ...
    'topo_vector').topo_vector;

T = 11; % 10 snapshots
interval = 1; % interval : 0.2s
unit = 0.1; % 0.1s
start1 = 77; % 46
start2 = 8;

fig = figure(1);
fig.Position = [50 100 1500 400]/fs_factor;
clf;
hold on;
for t = 1:T
    subplot(2,T,t);
    topo = nan(200);
    topo(topo_idx_comb) = topo_vector1(interval*(t-1)+1+start1,:);
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    if mod(t,2)
        text(0.5,-0.2, sprintf('%04.2fs', interval*(t-1)*unit), ...
            'HorizontalAlignment','center', 'Units','normalized', ...
            'FontSize',12);
    end
    if t==6
        text(0.5,-0.4, 'Time', ...
            'HorizontalAlignment','center', 'Units','normalized', ...
            'FontSize',16);
    end
end

for t = 1:T
    subplot(2,T,t+T);
    topo = nan(200);
    topo(topo_idx_comb) = topo_vector2(interval*(t-1)+1+start2,:);

    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    if mod(t,2)
        text(0.5,-0.2, sprintf('%04.2fs', interval*(t-1)*unit), ...
            'HorizontalAlignment','center', 'Units','normalized', ...
            'FontSize',12);
    end
    if t==6
        text(0.5,-0.4, 'Time', ...
            'HorizontalAlignment','center', 'Units','normalized', ...
            'FontSize',16);
    end
end
annotation('arrow', [0.12 0.92], [0.61 0.61], 'Units','normalized');
annotation('arrow', [0.12 0.92], [0.135 0.135], 'Units','normalized');

text(-14.3,2.45,sprintf('Conscious\n(resting)'), 'FontSize',14, ...
    'Units','normalized', 'Rotation',90, 'HorizontalAlignment','center');
text(-14.3,0.45,sprintf('Unconscious\n(suppression)'), 'FontSize',14, ...
    'Units','normalized', 'Rotation',90, 'HorizontalAlignment','center');

cb = colorbar;
cb.Location = 'north';
cb.Position = [0.758 0.43 0.1 0.025];
cb.Ticks = [-0.5 0.5];
cb.TickLabels = {'Lag', 'Lead'};
cb.FontSize = 10;
cb.Label.String = 'Relative phase';
cb.Label.FontSize = 12;

axes(fig, 'Position',[0.753 0.42 0.11 0.15], 'Units','normalized');
axis off;
patch([0 0 1 1], [0 1 1 0], 'w', 'FaceColor','none', ...
    'EdgeColor','k');
xlim([0 1]);
ylim([0 1]);

tic;
exportgraphics(fig, 'main_figures_UofM/png/snapshot_of_topoplot.png', 'Resolution',450);
exportgraphics(fig, 'main_figures_UofM/eps/snapshot_of_topoplot.eps', 'Resolution',450);
savefig(fig, 'main_figures_UofM/fig/snapshot_of_topoplot.fig');
fig.Color = 'w';
export_fig(fig, 'main_figures_UofM/pdf/snapshot_of_topoplot.pdf', '-pdf', '-opengl', '-r450');
toc;

%% bootstrap the UofM data
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

n_bootstrp = 1000;

load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/clustering_w_1min/properties_all.mat'],dropbox_path))
[s1, s2, s3, s4] = size(occur); n_pnt = s2*s3;

rng(240162,'twister'); % random number seed : my SKKU ID

occur_bootstrp   = nan([s1,n_pnt,n_bootstrp,s4]);
occur_w_bootstrp = nan([s1,n_pnt,n_bootstrp,s4]);
dwell_bootstrp   = nan([s1,n_pnt,n_bootstrp,s4]);
dwell_w_bootstrp = nan([s1,n_pnt,n_bootstrp,s4]);
for k = 1:s1
    for st = 1:n_state
        occur_tmp   = occur  (k,:,:,st);
        occur_w_tmp = occur_w(k,:,:,st);
        dwell_tmp   = dwell  (k,:,:,st);
        dwell_w_tmp = dwell_w(k,:,:,st);
        val_idx = find(~isnan(occur_tmp));
        n_idx = length(val_idx);
        [~, occur_btstrp_idx  ] = bootstrp(n_bootstrp,@mean,occur_tmp  (val_idx));
        [~, occur_w_btstrp_idx] = bootstrp(n_bootstrp,@mean,occur_w_tmp(val_idx));
        [~, dwell_btstrp_idx  ] = bootstrp(n_bootstrp,@mean,dwell_tmp  (val_idx));
        [~, dwell_w_btstrp_idx] = bootstrp(n_bootstrp,@mean,dwell_w_tmp(val_idx));
        for nb = 1:n_bootstrp
            occur_bootstrp  (k,1:n_idx,nb,st) = occur_tmp  (val_idx(occur_btstrp_idx  (:,nb)));
            occur_w_bootstrp(k,1:n_idx,nb,st) = occur_w_tmp(val_idx(occur_w_btstrp_idx(:,nb)));
            dwell_bootstrp  (k,1:n_idx,nb,st) = dwell_tmp  (val_idx(dwell_btstrp_idx  (:,nb)));
            dwell_w_bootstrp(k,1:n_idx,nb,st) = dwell_w_tmp(val_idx(dwell_w_btstrp_idx(:,nb)));
        end
    end
end

desc_columns = {'K', 'n_pnt', 'n_bootstrp', 'n_state'};
save('references/properties_bootstrp_all.mat', ...
    'desc_columns', 'state_list', 'n_state', 'n_bootstrp', 'n_pnt', ...
    'occur_bootstrp', 'occur_w_bootstrp', 'dwell_bootstrp', 'dwell_w_bootstrp', '-v7.3');

%% probability density function over brain states (load & draw : PDF) violin plot
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

K = 4;
noto_thr = 0.15; noto_thr = strrep(num2str(noto_thr),'.','p');
y_label = 'Probability density (%)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0 63]; y_height = diff(y_range);

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

% load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
%     'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/mode_distribution.mat'],dropbox_path));
% [s1, s2, s3, s4] = size(data);
% % K, n_slice, n_sub, n_state
% data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;
% % K, n_point, n_state -> K, n_state, n_point
load('references/properties_bootstrp_all.mat');
data = mean(occur_bootstrp,2,'omitnan');
[s1, s2, s3, s4] = size(data);
% K, n_bootstrp, n_state -> K, n_state, n_bootstrp
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;

state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
state_label_list = {'Eyes-closed', 'LOC', 'Burst', 'Suppression', 'ROC'};
n_state = length(state_list);

for k = 1:K+1
    x_list = (1:n_state)+(k-1)*(n_state+1)-1;
    y_mean = mean(data(k,:,:),3,'omitnan');
    plot(x_list+1, mean(data(k,:,:),3,'omitnan'), ':', ...
        'LineWidth',4, 'Color',c_list(k,:));
    % data shape to use da--plot function : y-axes by x-axes
    % for example, if you want to draw the mode PDF across level of
    % consciousness, the shape of input should be that 1st dimension of matrix
    % indicates data points and 2nd dimension shows brain states (90 by 5)
%     % box plot
%     h = daboxplot(squeeze(data(k,:,:))', 'fill',1, 'colors',c_list(k,:), ...
%         'mean',1, 'boxalpha',0.7, 'boxwidth',2, 'jitter',0);
    % violin plot
    h = daviolinplot(squeeze(data(k,:,:))', 'violin','full', 'colors',c_list(k,:), ...
        'violinalpha',0.7, 'boxwidth',1.6, 'outsymbol','rx');
    for st = 1:n_state
        h.bx(st).XData = h.bx(st).XData+(k-1)*(n_state+1);
        h.md(st).XData = h.md(st).XData+(k-1)*(n_state+1);
%         h.mn(st).XData = h.mn(st).XData+(k-1)*(n_state+1); % when drawing as box plot
        h.ds(st).XData = h.ds(st).XData+(k-1)*(n_state+1); % when drawing as violin plot
        h.ot(st).XData = h.ot(st).XData+(k-1)*(n_state+1);
        for j = 1:size(h.wh,3); h.wh(st,1,j).XData = h.wh(st,1,j).XData+(k-1)*(n_state+1); end
    end
    for st = 1:n_state; h.mn(st).YData = [y_mean(st) y_mean(st)]; end
end

load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/clustering_w_1min/results_100ms_bootstrp/occur.mat'],dropbox_path), 'pairs', 'pval');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(pairs,1)
        if pval{k}(te) <= 0.05
            x_l = pairs(te,1); x_r = pairs(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if pval{k}(te) <= 0.001; ast = ast_arr{3};
            elseif pval{k}(te) <= 0.01; ast = ast_arr{2};
            else; ast = ast_arr{1};
            end
            text((x_l+x_r)*0.5+(k-1)*(n_state+1), y_range(1)+y_height*0.575+ast_y_arr(te)*y_height*0.03, ast, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end

set(gca, 'FontSize',14);
xlabel('State', 'FontSize',18);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list state_label_list state_label_list]);
xlim([0 (n_state+1)*(K+1)]);
ylim(y_range);
ylabel(y_label, 'FontSize',18);

ah = axes(fig, 'Position',[0.35 0.626 0.3 0.15]);
patch([0.6 5.4 5.4 4 0.6], ...
    [0 0 y_height 0 y_height]/10, 'k', 'FaceAlpha',0.7, 'LineStyle','none');
set(ah, 'FontSize',12)
xticks(1:n_state);
xticklabels(state_label_list);
xlim([0.2 n_state+0.8])
yticks([]);
ylim([0 10]);
box on;
text(3,7,'Level of consciousness', 'FontSize',15, 'FontWeight','bold', ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center');

tic;
exportgraphics(fig, sprintf('main_figures_UofM/png/mode_distribution_w_noto_%s_violin.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('main_figures_UofM/eps/mode_distribution_w_noto_%s_violin.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('main_figures_UofM/fig/mode_distribution_w_noto_%s_violin.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('main_figures_UofM/pdf/mode_distribution_w_noto_%s_violin.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : dwell time) violin plot
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

K = 4;
noto_thr = 0.15; noto_thr = strrep(num2str(noto_thr),'.','p');
y_label = 'Dwell time (s)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0.08 0.32]; y_height = diff(y_range);

% to draw topographic map of K+1 centroids
load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', 'topo_idx_comb');
C_comb.EOEC = permute(C_comb.EOEC,[2,3,1]);
load('references/topoplot_parameters.mat');

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

% load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
%     'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/dwell_time.mat'],dropbox_path));
% [s1, s2, s3, s4] = size(data);
% % K, n_slice, n_sub, n_state
% data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);
% % K, n_point, n_state
load('references/properties_bootstrp_all.mat');
data = mean(dwell_bootstrp,2,'omitnan');
[s1, s2, s3, s4] = size(data);
% K, n_bootstrp, n_state -> K, n_state, n_bootstrp
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);

state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
state_label_list = {'Eyes-closed', 'LOC', 'Burst', 'Suppression', 'ROC'};
n_state = length(state_list);

for k = 1:K+1
    x_list = (1:n_state)+(k-1)*(n_state+1)-1;
    y_mean = mean(data(k,:,:),3,'omitnan');
    plot(x_list+1, mean(data(k,:,:),3,'omitnan'), ':', ...
        'LineWidth',4, 'Color',c_list(k,:));
    % data shape to use da--plot function : y-axes by x-axes
    % for example, if you want to draw the mode PDF across level of
    % consciousness, the shape of input should be that 1st dimension of matrix
    % indicates data points and 2nd dimension shows brain states (90 by 5)
%     % box plot
%     h = daboxplot(squeeze(data(k,:,:))', 'fill',1, 'colors',c_list(k,:), ...
%         'mean',1, 'boxalpha',0.7, 'boxwidth',2, 'jitter',0);
    % violin plot
    h = daviolinplot(squeeze(data(k,:,:))', 'violin','full', 'colors',c_list(k,:), ...
        'violinalpha',0.7, 'boxwidth',1.6, 'outsymbol','rx');
    for st = 1:n_state
        h.bx(st).XData = h.bx(st).XData+(k-1)*(n_state+1);
        h.md(st).XData = h.md(st).XData+(k-1)*(n_state+1);
%         h.mn(st).XData = h.mn(st).XData+(k-1)*(n_state+1); % when drawing as box plot
        h.ds(st).XData = h.ds(st).XData+(k-1)*(n_state+1); % when drawing as violin plot
        h.ot(st).XData = h.ot(st).XData+(k-1)*(n_state+1);
        for j = 1:size(h.wh,3); h.wh(st,1,j).XData = h.wh(st,1,j).XData+(k-1)*(n_state+1); end
    end
    for st = 1:n_state; h.mn(st).YData = [y_mean(st) y_mean(st)]; end
%     pl(k) = h.bx(1); % when drawing as box plot
    pl(k) = h.ds(1);  % when drawing as violin plot
end

load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/clustering_w_1min/results_100ms_bootstrp/dwell.mat'],dropbox_path), 'pairs', 'pval');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(pairs,1)
        if pval{k}(te) <= 0.05
            x_l = pairs(te,1); x_r = pairs(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.66+ast_y_arr(te)*y_height*0.03, 'k-');
            if pval{k}(te) <= 0.001; ast = ast_arr{3};
            elseif pval{k}(te) <= 0.01; ast = ast_arr{2};
            else; ast = ast_arr{1};
            end
            text((x_l+x_r)*0.5+(k-1)*(n_state+1), y_range(1)+y_height*0.665+ast_y_arr(te)*y_height*0.03, ast, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end

set(gca, 'FontSize',14);
xlabel('State', 'FontSize',18);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list state_label_list state_label_list]);
xlim([0 (n_state+1)*(K+1)]);
ylim(y_range);
ylabel(y_label, 'FontSize',18);

gap_margin = 0.184;
legend_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
for k = 1:K
    ah = axes('position',get(gca,'position'),'visible','off');
    legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
        'Location',[0.107+(k-1)*gap_margin 0.864 0.07 0.07], 'Units','normalized');
    
    ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.64 0.08 0.22], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end

k = K+1;
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah, pl(k), legend_label{k}, 'FontSize',12, ...
    'Location',[0.107+(k-1)*gap_margin 0.825 0.07 0.07], 'Units','normalized');

ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.6 0.08 0.22], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);

tic;
exportgraphics(fig, sprintf('main_figures_UofM/png/mode_dwell_time_w_noto_%s_violin.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('main_figures_UofM/eps/mode_dwell_time_w_noto_%s_violin.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('main_figures_UofM/fig/mode_dwell_time_w_noto_%s_violin.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('main_figures_UofM/pdf/mode_dwell_time_w_noto_%s_violin.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% inverse participation ratio (save)
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

tic;
load('../../../UM_data/pca/pca_st_EC_all.mat', 'coeff');
toc;
load('references/UM_info_ch128.mat');
load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/topo_vector'], dropbox_path);
state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
state_label_list = {'Eyes Closed', 'LOC', 'Busrt', 'Suppression', 'ROC'};
n_state = length(state_list);
n_slice = 5; % 5min -> 1min * 5
resol = 100; % 100ms
len_slice = 60/(resol/1e3);
K = 4;

iPR_all = nan([n_state,n_sub,n_slice]);
explained_all = nan([31739,n_state,n_sub,n_slice]);
score_all = cell([n_state,1]);
IDX_all = cell([n_state,1]);
IDX_all_w_noto = cell([n_state,1]);

load('references/mode_dist_profiles_w_noto_0p15_250306.mat', 'thr_arr');

for sub = 1:n_sub
    for st = 1:n_state
        state = state_list{st};
        try
            disp(sprintf('%s - %s', subjname{sub},state))
            load(sprintf('%s/topo_vector_st_%s_%s.mat',load_path,state,subjname{sub}));
            for sl = 1:n_slice
                idx_list = 1+(sl-1)*len_slice:sl*len_slice;
                score = topo_vector(idx_list,:)*coeff;
                explained = var(score,[],1)';
                explained = explained/sum(explained)*100;
                % iPR = \sum{\lambda^2} / (\sum{\lambda})^2
                iPR_all(st,sub,sl) = sum(explained'*explained)/sum(explained)^2;
                explained_all(:,st,sub,sl) = explained;
                if (st==1) || (st==4)
                    score_all{st} = cat(1,score_all{st},score(:,1:20));
                end
            end

            if (st==1) || (st==4)            
                load(sprintf('%s/../clustering/regr_st_%s_%s_n_nulls_1000.mat', ...
                    load_path,state,subjname{sub}), 'IDX');
                load(sprintf('%s/../regression_mae/regr_st_%s_%s.mat', ...
                    load_path,state,subjname{sub}), 'epsilon');
                IDX_all{st} = cat(1,IDX_all{st},IDX);
                IDX(epsilon>thr_arr(sub)) = K+1;
                IDX_all_w_noto{st} = cat(1,IDX_all_w_noto{st},IDX);
            end
        end
    end
end

tic;
save('references/iPR_matrix.mat', ...
    'state_list', 'state_label_list', 'n_sub', 'n_state', ...
    'explained_all', 'iPR_all', '-v7.3');
save(['references/' ...
    'score_all_matrix.mat'], ...
    'state_list', 'state_label_list', 'n_sub', 'n_state', ...
    'score_all', 'IDX_all', 'IDX_all_w_noto', '-v7.3');
toc;
% 
% [s1, s2, s3] = size(iPR_all);
% iPR_all_v2 = reshape(iPR_all, [s1, s2*s3]);

%% bootstrap the UofM data iPR
clear, close all;

n_bootstrp = 1000;

load('references/iPR_matrix.mat', ...
    'state_list', 'n_state', 'iPR_all');
% n_state, n_sub, n_slice -> n_slice, n_sub, n_state
iPR_all = permute(iPR_all,[3,2,1]);
[s1, s2, s3] = size(iPR_all); n_pnt = s1*s2;

rng(240162+1,'twister'); % random number seed : my SKKU ID +1 (to make a difference from before)

iPR_all_bootstrp   = nan([n_pnt,n_bootstrp,s3]);
for st = 1:n_state
    iPR_tmp = iPR_all(:,:,st);
    val_idx = find(~isnan(iPR_tmp));
    n_idx = length(val_idx);
    [~, bootstrp_idx  ] = bootstrp(n_bootstrp,@mean,iPR_tmp(val_idx));
    for nb = 1:n_bootstrp
        iPR_all_bootstrp(1:n_idx,nb,st) = iPR_tmp(val_idx(bootstrp_idx(:,nb)));
    end
end

desc_columns = {'n_pnt', 'n_bootstrp', 'n_state'};
save('references/iPR_matrix_bootstrp_all.mat', ...
    'desc_columns', 'state_list', 'n_state', 'n_bootstrp', 'n_pnt', ...
    'iPR_all_bootstrp', '-v7.3');


%% inverse participation ratio (draw)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

y_range = [0.09 0.46]; y_height = diff(y_range);

load('references/iPR_matrix.mat', 'explained_all', 'state_label_list');
% 31739,n_state,n_sub,n_slice -> 31739, n_state, n_pnt
[s1, s2, s3, s4] = size(explained_all);
explained_all = reshape(explained_all, [s1,s2,s3*s4]);
load('references/iPR_matrix_bootstrp_all.mat');
iPR_all_bootstrp = permute(iPR_all_bootstrp,[3,2,1]);
desc_columns = desc_columns([3,2,1]);
% load_path = sprintf('%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/band_[8-12]/100ms/clustering/anova',dropbox_path);
% load(sprintf('%s/inverse_participation_ratio.mat',load_path), 'results');
ast_arr = {'*', '**', '***'};
load('references/score_all_matrix.mat', 'score_all', 'IDX_all');
cum_explained_all = cumsum(explained_all,1);

fig = figure(1);
fig.Position = [100 300 1500 400]/fs_factor;
clf;
subplot(6,3,[3 6 9 12 15]);
hold on; grid on;

iPR_violin = mean(iPR_all_bootstrp,3,'omitnan');
plot(1:n_state, mean(iPR_violin,2,'omitnan'), ':', ...
    'LineWidth',4, 'Color',hex2rgb('#A2142F'));
daviolinplot(iPR_violin', 'violin','full', 'colors',hex2rgb('#A2142F'), ...
        'violinalpha',0.7, 'boxwidth',1.8, 'outsymbol','rx');

load(sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/clustering_w_1min/results_100ms_bootstrp/iPR.mat'],dropbox_path), 'pairs', 'pval');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for te = 1:size(pairs,1)
    if pval(te) <= 0.05
        x_l = pairs(te,1); x_r = pairs(te,2); 
        plot([x_l+0.1 x_r-0.1], y_range(1)+[y_height y_height].*0.85+ast_y_arr(te)*y_height*0.03, 'k-');
        if pval(te) <= 0.001; ast = ast_arr{3};
        elseif pval(te) <= 0.01; ast = ast_arr{2};
        else; ast = ast_arr{1};
        end
        text((x_l+x_r)*0.5, y_range(1)+y_height*0.855+ast_y_arr(te)*y_height*0.03, ast, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end
xticks(1:n_state);
xticklabels(state_label_list);
xlabel('State', 'FontSize',14);
ylabel('inverse Participation Ratio', 'FontSize',14);
xlim([0 n_state+1]);
ylim(y_range);

n_ax = 20;
% ========================== EC explained ===========================
subplot(2,3,1);
hold on; grid on;

y = mean(cum_explained_all(:,1,:),3,'omitnan');
y = y(1:n_ax);
se_y = std(cum_explained_all(:,1,:),0,3,'omitnan')./sqrt(sum(~isnan(cum_explained_all(1,1,:)),3));
se_y = se_y(1:n_ax);
patch([1:n_ax flip(1:n_ax)], [y-se_y; flip(y+se_y)], ...
    'b', 'lineStyle','none', 'Facealpha',0.2);
plot(y, '.-');
yline(75, 'k--');
xlim([1 n_ax]);
ylim([0 100]);
xlabel('# of principal components', 'FontSize',12);
ylabel('Explainability (%)', 'FontSize',12);
legend({'standard error', 'average', '75% guideline'}, ...
    'Location','southeast', 'FontSize',12);

% ========================== Suppr. explained ===========================
subplot(2,3,4);
hold on; grid on;

y = mean(cum_explained_all(:,4,:),3,'omitnan');
y = y(1:n_ax);
se_y = std(cum_explained_all(:,4,:),0,3,'omitnan')./sqrt(sum(~isnan(cum_explained_all(1,4,:)),3));
se_y = se_y(1:n_ax);
patch([1:n_ax flip(1:n_ax)], [y-se_y; flip(y+se_y)], ...
    'b', 'lineStyle','none', 'Facealpha',0.2);
plot(y, '.-');
yline(75, 'k--');
xlim([1 n_ax]);
ylim([0 100]);
xlabel('# of principal components', 'FontSize',12);
ylabel('Explainability (%)', 'FontSize',12);
legend({'standard error', 'average', '75% guideline'}, ...
    'Location','southeast', 'FontSize',12);

K = 4;
c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;

% ========================== EC scatter ===========================
scatter_legend_list = {'Mode1', 'Mode2', 'Mode3', 'Mode4'};
subplot(2,3,2);
hold on; grid on;

st = 1;
for k = 1:K
    idx_tmp = (IDX_all{st}==k);
    sc(k) = scatter(score_all{st}(idx_tmp,2),score_all{st}(idx_tmp,1), 15, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor',c_list(k,:), ...
        'MarkerEdgeAlpha',0.6, 'MarkerFaceAlpha',0.6);
end
xlim([-200 200]);
ylim([-200 200]);
xline(0, 'LineWidth',2);
yline(0, 'LineWidth',2);
xlabel('2^{nd} PC', 'FontSize',12);
ylabel('1^{st} PC', 'FontSize',12);
legend(sc, scatter_legend_list, 'Location','northeastoutside', ...
    'FontSize',12);
text(-0.3, 1.14, 'Eyes-closed state', 'FontSize',15, 'FontWeight','normal', ...
    'Units','normalized', 'HorizontalAlignment','center');

% ========================== Suppr. scatter ===========================
subplot(2,3,5);
hold on; grid on;

st = 4;
for k = 1:K
    idx_tmp = (IDX_all{st}==k);
    sc(k) = scatter(score_all{st}(idx_tmp,2),score_all{st}(idx_tmp,1), 15, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor',c_list(k,:), ...
        'MarkerEdgeAlpha',0.6, 'MarkerFaceAlpha',0.6);
end
xlim([-200 200]);
ylim([-200 200]);
xline(0, 'LineWidth',2);
yline(0, 'LineWidth',2);
xlabel('2^{nd} PC', 'FontSize',12);
ylabel('1^{st} PC', 'FontSize',12);
legend(sc, scatter_legend_list, 'Location','northeastoutside', ...
    'FontSize',12);
text(-0.3, 1.14, 'Suppression state', 'FontSize',15, 'FontWeight','normal', ...
    'Units','normalized', 'HorizontalAlignment','center');

tic;
exportgraphics(fig, 'main_figures_UofM/png/PCA_result_w_iPR.png', 'Resolution',450);
exportgraphics(fig, 'main_figures_UofM/eps/PCA_result_w_iPR.eps', 'Resolution',450);
savefig(fig, 'main_figures_UofM/fig/PCA_result_w_iPR.fig');
fig.Color = 'w';
export_fig(fig, 'main_figures_UofM/pdf/PCA_result_w_iPR.pdf', '-pdf', '-opengl', '-r450');
toc;

%% averaged topomap of each mode with noto (save)
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

load('references/UM_info_ch128.mat');
load('references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'topo_idx_UM2comb','topo_idx_comb');
noto_thr = 0.15; noto_thr = strrep(num2str(noto_thr),'.','p');
load(sprintf('references/mode_dist_profiles_w_noto_%s_250306.mat',noto_thr), ...
    'thr_arr','thr_ratio');
load_path = sprintf('%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/band_[8-12]/100ms/',dropbox_path);
state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
% state = 'EC';
% state = 'DS';

K = 4;
unit = 0.1;

for st = [1 4]
    tic;
    state = state_list{st};

    topo_vector_all = zeros([length(topo_idx_comb), K+1]);
    IDX_all = zeros([1,K+1]);
    for idx = 1:n_sub    
        try
            load([load_path '/topo_vector/topo_vector_st_' ...
                state '_' subjname{idx} '.mat']);
            load([load_path '/regression/regr_st_' ...
                state '_' subjname{idx} '.mat'], 'beta', 'epsilon');
            IDX = cal_regression_clustering(beta,K,unit,1);
            IDX(mean(abs(epsilon),1)>thr_arr(idx)) = 5;
    
            IDX_all = IDX_all + histcounts(IDX,0.5:K+1.5,'Normalization','count');
        
            for k = 1:K+1
                topo_vector_all(:,k) = topo_vector_all(:,k) + ...
                    sum(topo_vector(IDX==k,:),1)';
            end
        end
    end
    topo_vector_all = topo_vector_all./IDX_all;
    
    save(['references/mean_topo_vector_for_each_k_in_st_' state '_w_noto_' noto_thr '.mat'], ...
        'K', 'topo_vector_all', 'IDX_all', 'unit', ...
        'thr_arr', 'thr_ratio', '-v7.3');
    toc;
end

%% draw transtion matrix (version 2)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

load_path = sprintf('%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/band_[8-12]/clustering',dropbox_path);
load('../../../results_yjp/UM_info_ch128.mat');

load('references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'topo_idx_comb');
load('references/topoplot_parameters.mat');

noto_thr = 0.15;
% for noto_thr = 0.05:0.05:0.3
    noto_thr_s = strrep(num2str(noto_thr),'.','p');
    load(sprintf('references/mode_dist_profiles_w_noto_%s_250306.mat',noto_thr_s), ...
        'trans_dist_all', 'readme_last_column_trans', 'readme_axis');
    
    K = 4; unit = 0.1;
    
    % load transition matrix
    % K+1, K+1, n_sub, n_slice, state
    trans_mat_all = permute(trans_dist_all(:,:,[1 4],:,:,1),[1,2,4,5,3]);

    
    % calculate the differences and p-values
    % (DS - EC)/EC * 100
    trans_mat_mean_diff = ...
        (mean(trans_mat_all(:,:,:,:,2),[3,4],'omitnan') ...
        - mean(trans_mat_all(:,:,:,:,1),[3,4],'omitnan')) ...
        ./ mean(trans_mat_all(:,:,:,:,1),[3,4],'omitnan') * 100;
    
    ttest2_mat = nan(K+1);
    for i = 1:K+1
        for j = 1:K+1
            val1 = squeeze(trans_mat_all(i,j,:,:,1)); val1 = val1(~isnan(val1));
            val2 = squeeze(trans_mat_all(i,j,:,:,2)); val2 = val2(~isnan(val2));
            [~,ttest2_mat(i,j)] = ttest2(val1,val2);
        end
    end
    
    
    margin_height = 0.16;
    margin_width = 0.04;
    
    fig = figure(1);
    fig.Position = [100 100 1600 500]/fs_factor;
    clf;
    
    % draw heatmap
    ax1 = subplot_tight(1,3,1,[margin_height,margin_width]);
    hold on; grid on;
    imagesc(-log10(ttest2_mat)', [1,3]);
    set(gca, 'YDir','reverse');
%     colormap(gca, flip(slanCM('bone')));
    cb = colorbar;
    cb.Label.String = 'p-value';
    cb.Label.FontSize = 18;
    cb.Ticks = 1:3;
    cb.TickLabels = 10.^-(1:3);
    set(gca, 'FontSize',15);
    xlabel('Mode (to)', 'FontSize',18);
    xticks(1:K+1); xlim([0.5,K+1.5]);
    ylabel('Mode (from)', 'FontSize',18);
    yticks(1:K+1); ylim([0.5,K+1.5]);
    [xTxt, yTxt] = ndgrid(1:5, 1:5); 
    
    flatten_trans_mat_mean_diff = reshape(trans_mat_mean_diff,1,[]);
    txt_cell = cell([1,length(flatten_trans_mat_mean_diff)]);
    for tx = 1:length(txt_cell); txt_cell{tx} = sprintf('%.2f',flatten_trans_mat_mean_diff(tx)); end
    
    text(xTxt(:), yTxt(:), txt_cell, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
        'FontSize',15, 'Color','w');
    xline(0.5:K+1.5); yline(0.5:K+1.5);
    title('Relative Change (%)', 'FontSize',16, 'FontWeight','normal');
    
    % ---------------------------------------------------------------
    state = 'EC';
    trans_mean = squeeze(mean(trans_mat_all(:,:,:,:,1),[3,4],'omitnan'));
    G = digraph(trans_mean);
    
    ax2 = subplot_tight(1,3,2,[margin_height,margin_width]);
    ax2.Position(1) = ax2.Position(1)+0.04;
    ax2.Position(2) = ax2.Position(2)-0.1;
    ax2.Position(4) = ax2.Position(4)+0.08;
    axis off;
    
    x = [0 -1 1  0 1.1]; y = [1  0 0 -1 -1.1];
    load(['references/mean_topo_vector_for_each_k_in_st_' state '_w_noto_' noto_thr_s '.mat']);
    x_n = [0.114 0.027 0.202 0.114 0.210]; 
    y_n = [0.526 0.296 0.296 0.065 0.046];
    xwidth = 0.052; ywidth = 0.19;
    % 1:[0.114 0.526], 2:[0.027 0.296], 3:[0.202 0.296]
    % 4:[0.114 0.065], 5:[0.210 0.046]
    for k = 1:K+1
        ax_in = axes('Position', ...
            [ax2.Position(1:2)+[x_n(k) y_n(k)] xwidth ywidth], ...
            'Units','normalized');
        hold on;
    
        topo = nan(200);
        topo(topo_idx_comb) = topo_vector_all(:,k);
        
        topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
            'scatter',0,'line',0);
        colormap(ax_in, slanCM('jet'));
        clim(ax_in, [-0.75 0.75]);
    end
    
    ax2 = subplot_tight(1,3,2,[margin_height,margin_width]);
    ax2.Position(1) = ax2.Position(1)+0.04;
    ax2.Position(2) = ax2.Position(2)-0.1;
    ax2.Position(4) = ax2.Position(4)+0.08;
    hold on; grid on;

    xlim([-1.6 1.6]); ylim([-1.7 1.6]);
    axis off;

    p = plot(ax2, G);
    p.XData = x; p.YData = y;
    
    p.Marker = 'none';
    p.NodeLabel = repmat({''},1,numnodes(G));
    p.LineWidth = 10*G.Edges.Weight;
    p.ArrowSize = 0*G.Edges.Weight;
    p.EdgeCData = G.Edges.Weight; p.EdgeAlpha = 1;
    
    G = digraph(trans_mean-diag(diag(trans_mean)));
    p = plot(ax2, G);
    p.XData = x; p.YData = y;
    
    p.Marker = '.'; 
    p.MarkerSize = 25*ones([1,K+1]); p.NodeColor = [60 60 70]./255;
    p.NodeLabel = repmat({''},1,numnodes(G));
    p.LineWidth = 10*G.Edges.Weight;
    p.ArrowSize = 60*G.Edges.Weight; p.ArrowPosition = 0.65;
%     colormap(ax2, flipud(gray)); clim(ax2, [0,0.4]);
    p.EdgeCData = G.Edges.Weight; p.EdgeAlpha = 1;
    
    text(p.XData,p.YData+0.47,{'Mode 1','Mode 2','Mode 3','Mode 4','Other'}, ...
        'FontSize',14, 'HorizontalAlignment','center', 'Color',[0 0 0]);
    text(0, ax2.YLim(2)+0.2, 'Eyes-closed state', ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
        'FontSize',16, 'FontWeight','normal');
    
    % ---------------------------------------------------------------
    state = 'DS';
    trans_mean = squeeze(mean(trans_mat_all(:,:,:,:,2),[3,4],'omitnan'));
    G = digraph(trans_mean);
    
    ax3 = subplot_tight(1,3,3,[margin_height,margin_width]);
    ax3.Position(2) = ax3.Position(2)-0.1;
    ax3.Position(4) = ax3.Position(4)+0.08;
    axis off;
    
    x = [0 -1 1  0 1.1]; y = [1  0 0 -1 -1.1];
    load(['references/mean_topo_vector_for_each_k_in_st_' state '_w_noto_' noto_thr_s '.mat']);
    x_n = [0.114 0.027 0.202 0.114 0.210]; 
    y_n = [0.526 0.296 0.296 0.065 0.046];
    xwidth = 0.052; ywidth = 0.19;
    % 1:[0.114 0.526], 2:[0.027 0.296], 3:[0.202 0.296]
    % 4:[0.114 0.065], 5:[0.210 0.046]
    for k = 1:K+1
        ax_in = axes('Position', ...
            [ax3.Position(1:2)+[x_n(k) y_n(k)] xwidth ywidth], ...
            'Units','normalized');
        hold on;
    
        topo = nan(200);
        topo(topo_idx_comb) = topo_vector_all(:,k);
        
        topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
            'scatter',0,'line',0);
        colormap(ax_in, slanCM('jet'));
        clim(ax_in, [-0.75 0.75]);
    end

    ax3 = subplot_tight(1,3,3,[margin_height,margin_width]);
    ax3.Position(2) = ax3.Position(2)-0.1;
    ax3.Position(4) = ax3.Position(4)+0.08;
    hold on; grid on;

    xlim([-1.6 1.6]); ylim([-1.7 1.6]);
    axis off;
    
    p = plot(ax3, G);
    p.XData = x; p.YData = y;
    
    p.Marker = 'none';
    p.NodeLabel = repmat({''},1,numnodes(G));
    p.LineWidth = 10*G.Edges.Weight;
    p.ArrowSize = 0*G.Edges.Weight;
    p.EdgeCData = G.Edges.Weight; p.EdgeAlpha = 1;
    
    G = digraph(trans_mean-diag(diag(trans_mean)));
    p = plot(ax3, G);
    p.XData = x; p.YData = y;
    
    p.Marker = '.'; 
    p.MarkerSize = 25*ones([1,K+1]); p.NodeColor = [60 60 70]./255;
    p.NodeLabel = repmat({''},1,numnodes(G));
    p.LineWidth = 10*G.Edges.Weight;
    p.ArrowSize = 60*G.Edges.Weight; p.ArrowPosition = 0.65;
%     colormap(ax3, flipud(gray)); clim(ax3, [0,0.4]);
    p.EdgeCData = G.Edges.Weight; p.EdgeAlpha = 1;
    
    text(p.XData,p.YData+0.47,{'Mode 1','Mode 2','Mode 3','Mode 4','Other'}, ...
        'FontSize',14, 'HorizontalAlignment','center', 'Color',[0 0 0]);
    text(0, ax3.YLim(2)+0.2, 'Suppression state', ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
        'FontSize',16, 'FontWeight','normal');
    
    colormap(ax1, flip(slanCM('bone')));
    colormap(ax2, flipud(gray)); clim(ax2, [0,0.4]);
    colormap(ax3, flipud(gray)); clim(ax3, [0,0.4]);

    if noto_thr == 0.15
        tic;
        exportgraphics(fig, sprintf('main_figures_UofM/png/transition_matrix_w_noto_%s_v2.png',noto_thr_s), 'Resolution',450);
        exportgraphics(fig, sprintf('main_figures_UofM/eps/transition_matrix_w_noto_%s_v2.eps',noto_thr_s), 'Resolution',450);
        savefig(fig, sprintf('main_figures_UofM/fig/transition_matrix_w_noto_%s_v2.fig',noto_thr_s));
        fig.Color = 'w';
        export_fig(fig, sprintf('main_figures_UofM/pdf/transition_matrix_w_noto_%s_v2.pdf',noto_thr_s), '-pdf', '-opengl', '-r450');
        toc;
    else
        tic;
        exportgraphics(fig, sprintf('si_figures/png/transition_matrix_w_noto_%s_v2.png',noto_thr_s), 'Resolution',450);
        exportgraphics(fig, sprintf('si_figures/eps/transition_matrix_w_noto_%s_v2.eps',noto_thr_s), 'Resolution',450);
        savefig(fig, sprintf('si_figures/fig/transition_matrix_w_noto_%s_v2.fig',noto_thr_s));
        fig.Color = 'w';
        export_fig(fig, sprintf('si_figures/pdf/transition_matrix_w_noto_%s_v2.pdf',noto_thr_s), '-pdf', '-opengl', '-r450');
        toc;
    end
% end