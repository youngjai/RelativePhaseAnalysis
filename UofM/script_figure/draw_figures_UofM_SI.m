%% to find optimal the number of clusters, K, in UofM (save)
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

load('references/UM_info_ch128.mat');
load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref/' ...
    'band_[8-12]/100ms/topo_vector'],dropbox_path);

state_list = {'EC', 'DS'};
% state = 'EC';
for st = 1:2
    state = state_list{st};
    K_list = 1:10;
    Mat_eva = nan([length(K_list), n_sub]);
    for idx = 1:n_sub
        try
            load(sprintf('%s/topo_vector_st_%s_%s.mat',load_path,state,subjname{idx}));
            tic;
            disp(idx);
            eva = evalclusters(topo_vector,'kmeans', ...
                'DaviesBouldin', 'KList',K_list);
            toc;
            Mat_eva(:,idx) = eva.CriterionValues;
        end
    end
    save(['references/SI/eval_kmean_' state '.mat'], ...
        'Mat_eva', 'K_list', '-v7.3');
end

%% to find optimal the number of clusters, K, in UofM (load, draw)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

fig = figure(1);
fig.Position = [100 100 1200 450]/fs_factor;
clf;
% ================ Eyes-closed ====================
subplot(1,2,1);
hold on; grid on;
state = 'EC';
load(['references/SI/eval_kmean_' state '.mat']);
mean_eval = mean(Mat_eva,2,'omitnan');
se_eval = std(Mat_eva,0,2,'omitnan')./sqrt(sum(~isnan(Mat_eva(end,:))));
bar(K_list([1:3 5:end]), mean_eval([1:3 5:end]), ...
    'LineWidth',1.2, 'FaceAlpha',0.8, 'FaceColor','#0072BD');
bar(K_list(4), mean_eval(4), ...
    'LineWidth',1.2, 'FaceAlpha',0.8, 'FaceColor','#EDB120');
errorbar(K_list, mean_eval, se_eval, 'k', ...
    'CapSize',20, 'LineWidth',2, 'LineStyle','none');
set(gca, 'FontSize',12);
xlabel('K', 'FontSize',16);
ylabel('DBI', 'FontSize',16);
xlim([1 11]);
ylim([1.05 1.75]);
text(-0.16,1.054,'(A) Davies-Bouldin index in Eyes-closed state', ...
    'FontWeight','bold', 'FontSize',14, 'Units','normalized');

% ================ Eyes-closed ====================
subplot(1,2,2);
hold on; grid on;
state = 'DS';
load(['references/SI/eval_kmean_' state '.mat']);
mean_eval = mean(Mat_eva,2,'omitnan');
se_eval = std(Mat_eva,0,2,'omitnan')./sqrt(sum(~isnan(Mat_eva(end,:))));
bar(K_list([1:3 5:end]), mean_eval([1:3 5:end]), ...
    'LineWidth',1.2, 'FaceAlpha',0.8, 'FaceColor','#0072BD');
bar(K_list(4), mean_eval(4), ...
    'LineWidth',1.2, 'FaceAlpha',0.8, 'FaceColor','#EDB120');
errorbar(K_list, mean_eval, se_eval, 'k', ...
    'CapSize',20, 'LineWidth',2, 'LineStyle','none');
set(gca, 'FontSize',12);
xlabel('K', 'FontSize',16);
ylabel('DBI', 'FontSize',16);
xlim([1 11]);
ylim([1.05 2.9]);
text(-0.16,1.054,'(B) Davies-Bouldin index in Suppression state', ...
    'FontWeight','bold', 'FontSize',14, 'Units','normalized');

tic;
exportgraphics(fig, 'si_figures/png/davies_bouldin_index.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/davies_bouldin_index.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/davies_bouldin_index.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/davies_bouldin_index.pdf', '-pdf', '-opengl', '-r450');
toc;

%% to find optimal time scale of the relative phase in UofM datset (save)
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

load('references/UM_info_ch128.mat');
load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/' ...
    'preprocessed_wo_reref/band_[0p1-100]/'],dropbox_path);
% state = 'DS';
state_list = {'EC', 'DS'};
for st = 1:2
    state = state_list{st};
    eeglab nogui;
    
    band_list = [[1 2]; [2 3]; [3 4]; ...
        [4 6]; [6 8]; [8 10]; [10 12]; ...
        [12 16]; [16 20]; [20 24]; [24 28]; [28 32]; ...
        [32 40]; [40 48]];
    band_name_list = {'[ 1- 2]'; '[ 2- 3]'; '[ 3- 4]'; ...
        '[ 4- 6]'; '[ 6- 8]'; '[ 8-10]'; '[10-12]'; ...
        '[12-16]'; '[16-20]'; '[20-24]'; '[24-28]'; '[28-32]'; ...
        '[32-40]'; '[40-48]'};
    n_band = size(band_list,1);
    
    time_resol_list = [2:2:200, 250, 334, 500, 1000, 2000]; % unit: ms
    time_resol_list = time_resol_list/2;
    n_trl = length(time_resol_list);
    
    Mat_shannon_std = nan([n_trl, n_band, n_sub, 2]);
    name_list = {'time_resol(ms)', 'band', 'subj', 'index'};
    name_list_of_4th_column = {'shannon_entropy', 'std'};
    
    for idx = 1:n_sub
        try
            load([load_path subjname{idx} '_wo_badchan_and_wo_reref.mat']);
            st = strcmp(state_list,state);
            data = data_state{st}';
        
            for bl = 1:n_band
                band = band_list(bl,:);
                EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
                EEG = eeg_checkset(EEG);
                
                EEG = pop_eegfiltnew(EEG, 'locutoff',band(1), ...
                    'hicutoff',band(2));
                H = EEG.data'; clear EEG;
                rel_p = cal_rel_phase_v3(H);
            
                N = 1000;
                n_ch = size(rel_p,2);
                n_bins = linspace(-1,1,N+1);
            
                for trl = 1:n_trl
                    tr = time_resol_list(trl);
                    rel_p_w_m = moving_time_window_v2(rel_p,tr,tr,@(x) x);
                    
                    p_ch = zeros([n_ch, N]);
                    for c = 1:n_ch
                        p_ch(c,:) = histcounts(rel_p_w_m(:,c),n_bins, 'Normalization','probability');
                    end
                    Mat_shannon_std(trl,bl,idx,1) = mean(sum(-p_ch.*log2(p_ch),2,'omitnan')/log2(N));
                    Mat_shannon_std(trl,bl,idx,2) = mean(std(rel_p_w_m,0,1,'omitnan'),'omitnan');
                end
            end
        end
    end
    time_resol_list = 2.*time_resol_list;
    save(['references/SI/time_resol_check_' state '.mat'], ...
        'Mat_shannon_std', 'name_list', 'name_list_of_4th_column', ...
        'band_list', 'band_name_list', 'time_resol_list', '-v7.3');
end

%% to find optimal time scale of the relative phase in UofM datset (load, draw)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

state = 'EC';
load(['references/SI/time_resol_check_' state '.mat']);

fig = figure(3);
fig.Position = [100 100 1500 600]/fs_factor;
clf;
subplot_tight(1,2,1,0.14);
hold on; grid on;
text(-0.35,1.08,'(A)', 'FontSize',20, 'Units','normalized');
contourf(time_resol_list, 1:14, mean(Mat_shannon_std(:,:,:,1),3,'omitnan')', ...
    [0.6:0.04:0.96 0.969:0.008:1]);
[C,h] = contourf(time_resol_list, 1:14, mean(Mat_shannon_std(:,:,:,1),3,'omitnan')', ...
    [0.62:0.04:0.96 0.965:0.008:1], 'ShowText',true, 'LabelSpacing', 400, ...
    'FaceColor','none', 'LineWidth',2);
clabel(C,h, 'FontSize',14);
set(gca, 'XScale','log');
set(gca, 'FontSize',16);
xlabel('Time (ms)', 'FontSize',20);
xlim([time_resol_list(1), time_resol_list(end)]);
xticks([10 100 1000]);
ylabel('Frequency band (Hz)', 'FontSize',20);
ylim([1,length(band_name_list)]);
yticks(1:length(band_name_list));
yticklabels(band_name_list);
colormap(gca, slanCM('Blues'));
clb = colorbar;
clb.Label.String = "Shannon's entropy";
clb.Label.FontSize = 22;
clim([0.89 0.99]);

subplot_tight(1,2,2,0.14);
hold on; grid on;
text(-0.35,1.08,'(B)', 'FontSize',20, 'Units','normalized');
contourf(time_resol_list, 1:14, mean(Mat_shannon_std(:,:,:,2),3,'omitnan')', ...
    [0.025:0.05:0.55 0.59:0.01:0.7]);
[C,h] = contourf(time_resol_list, 1:14, mean(Mat_shannon_std(:,:,:,2),3,'omitnan')', ...
    [0:0.05:0.55 0.585:0.01:0.7], 'ShowText',true, 'LabelSpacing', 400, ...
    'FaceColor','none', 'LineWidth',2);
clabel(C,h, 'FontSize',14);
set(gca, 'XScale','log');
set(gca, 'FontSize',16);
xlabel('Time (ms)', 'FontSize',20);
xlim([time_resol_list(1), time_resol_list(end)]);
xticks([10 100 1000]);
ylabel('Frequency band (Hz)', 'FontSize',20);
ylim([1,length(band_name_list)]);
yticks(1:length(band_name_list));
yticklabels(band_name_list);
colormap(gca, slanCM('Greens'));
clb = colorbar;
clb.Label.String = 'Standard deviation';
clb.Label.FontSize = 22;
clim([0.37 0.63]);

tic;
exportgraphics(fig, sprintf('si_figures/png/time_resol_check_%s.png',state), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/time_resol_check_%s.eps',state), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/time_resol_check_%s.fig',state));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/time_resol_check_%s.pdf',state), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : PDF) violin plot (alpha 100ms, unweighted)
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
y_label = 'Normalized occurrence (%)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0 70]; y_height = diff(y_range);

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

data = load('references/properties_all_alpha_w_noto_0p15_100ms.mat','occur').occur;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;
% K, n_point, n_state -> K, n_state, n_point
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/weighted_mode_distribution.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

ah = axes(fig, 'Position',[0.35 0.65 0.3 0.15]);
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
exportgraphics(fig, sprintf('si_figures/png/mode_distribution_w_noto_%s_violin_alpha_100ms_unweigthed.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_distribution_w_noto_%s_violin_alpha_100ms_unweigthed.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_distribution_w_noto_%s_violin_alpha_100ms_unweigthed.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_distribution_w_noto_%s_violin_alpha_100ms_unweigthed.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : dwell time) violin plot (alpha 100ms, unweighted)
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
y_label = 'Normalized dwell time (s)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0.02 0.36]; y_height = diff(y_range);

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

data = load('references/properties_all_alpha_w_noto_0p15_100ms.mat','dwell').dwell;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);
% K, n_point, n_state
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/weighted_dwell_time.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

gap_margin = 0.184;
legend_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
for k = 1:K
    ah = axes('position',get(gca,'position'),'visible','off');
    legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
        'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');
    
    ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end

k = K+1;
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
    'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');

ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);

tic;
exportgraphics(fig, sprintf('si_figures/png/mode_dwell_time_w_noto_%s_violin_alpha_100ms_unweigthed.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_dwell_time_w_noto_%s_violin_alpha_100ms_unweigthed.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_dwell_time_w_noto_%s_violin_alpha_100ms_unweigthed.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_dwell_time_w_noto_%s_violin_alpha_100ms_unweigthed.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;


%% probability density function over brain states (load & draw : PDF) violin plot (alpha 100ms, weighted)
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
y_label = 'Normalized occurrence (%)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0 70]; y_height = diff(y_range);

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

data = load('references/properties_all_alpha_w_noto_0p15_100ms.mat','occur_w').occur_w;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;
% K, n_point, n_state -> K, n_state, n_point
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/weighted_mode_distribution.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

ah = axes(fig, 'Position',[0.35 0.65 0.3 0.15]);
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
exportgraphics(fig, sprintf('si_figures/png/mode_distribution_w_noto_%s_violin_alpha_100ms_weigthed.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_distribution_w_noto_%s_violin_alpha_100ms_weigthed.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_distribution_w_noto_%s_violin_alpha_100ms_weigthed.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_distribution_w_noto_%s_violin_alpha_100ms_weigthed.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : dwell time) violin plot (alpha 100ms, weighted)
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
y_label = 'Normalized dwell time (s)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0.02 0.36]; y_height = diff(y_range);

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

data = load('references/properties_all_alpha_w_noto_0p15_100ms.mat','dwell_w').dwell_w;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);
% K, n_point, n_state
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/100ms/clustering_w_1min/stats_anova_1min/weighted_dwell_time.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

gap_margin = 0.184;
legend_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
for k = 1:K
    ah = axes('position',get(gca,'position'),'visible','off');
    legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
        'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');
    
    ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end

k = K+1;
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
    'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');

ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);

tic;
exportgraphics(fig, sprintf('si_figures/png/mode_dwell_time_w_noto_%s_violin_alpha_100ms_weigthed.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_dwell_time_w_noto_%s_violin_alpha_100ms_weigthed.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_dwell_time_w_noto_%s_violin_alpha_100ms_weigthed.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_dwell_time_w_noto_%s_violin_alpha_100ms_weigthed.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : PDF) violin plot (alpha 20ms, unweighted)
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
y_range = [0 70]; y_height = diff(y_range);

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

data = load('references/properties_all_alpha_w_noto_0p15_20ms.mat','occur').occur;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;
% K, n_point, n_state -> K, n_state, n_point
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/20ms/clustering_w_1min/stats_anova_1min/occur.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

ah = axes(fig, 'Position',[0.35 0.65 0.3 0.15]);
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
exportgraphics(fig, sprintf('si_figures/png/mode_distribution_w_noto_%s_violin_alpha_20ms_unweighted.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_distribution_w_noto_%s_violin_alpha_20ms_unweighted.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_distribution_w_noto_%s_violin_alpha_20ms_unweighted.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_distribution_w_noto_%s_violin_alpha_20ms_unweighted.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : dwell time) violin plot (alpha 20ms, unweighted)
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
y_range = [0.02 0.36]; y_height = diff(y_range);

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

data = load('references/properties_all_alpha_w_noto_0p15_20ms.mat','dwell').dwell;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);
% K, n_point, n_state
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/20ms/clustering_w_1min/stats_anova_1min/dwell.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

gap_margin = 0.184;
legend_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
for k = 1:K
    ah = axes('position',get(gca,'position'),'visible','off');
    legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
        'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');
    
    ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end

k = K+1;
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
    'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');

ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);

tic;
exportgraphics(fig, sprintf('si_figures/png/mode_dwell_time_w_noto_%s_violin_alpha_20ms_unweighted.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_dwell_time_w_noto_%s_violin_alpha_20ms_unweighted.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_dwell_time_w_noto_%s_violin_alpha_20ms_unweighted.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_dwell_time_w_noto_%s_violin_alpha_20ms_unweighted.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : PDF) violin plot (alpha 20ms, weighted)
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
y_label = 'Normalized occurrence (%)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0 70]; y_height = diff(y_range);

fig = figure(1);
fig.Position = [100 50 1800 500]/fs_factor;
clf;
ax_main = axes(fig, 'Position',[0.05 0.25 0.92 0.45]);
hold on; grid on;

data = load('references/properties_all_alpha_w_noto_0p15_20ms.mat','occur_w').occur_w;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2])*100;
% K, n_point, n_state -> K, n_state, n_point
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/20ms/clustering_w_1min/stats_anova_1min/occur_w.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

ah = axes(fig, 'Position',[0.35 0.65 0.3 0.15]);
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
exportgraphics(fig, sprintf('si_figures/png/mode_distribution_w_noto_%s_violin_alpha_20ms_weighted.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_distribution_w_noto_%s_violin_alpha_20ms_weighted.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_distribution_w_noto_%s_violin_alpha_20ms_weighted.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_distribution_w_noto_%s_violin_alpha_20ms_weighted.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% probability density function over brain states (load & draw : dwell time) violin plot (alpha 20ms, weighted)
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
y_label = 'Normalized dwell time (s)';

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
y_range = [0.02 0.36]; y_height = diff(y_range);

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

data = load('references/properties_all_alpha_w_noto_0p15_20ms.mat','dwell_w').dwell_w;
[s1, s2, s3, s4] = size(data);
% K, n_slice, n_sub, n_state
data = reshape(data, [s1, s2*s3, s4]); data = permute(data,[1,3,2]);
% K, n_point, n_state
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
        'violinalpha',0.7, 'boxwidth',1.2, 'outsymbol','rx');
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
    'band_[8-12]/20ms/clustering_w_1min/stats_anova_1min/dwell_w.mat'],dropbox_path), 'results');
ast_arr = {'*', '**', '***'};
ast_y_arr = [0 1 3 5 0 2 4 0 1 0];
for k = 1:K+1
    for te = 1:size(results{k},1)
        if results{k}(te,6) <= 0.05
            x_l = results{k}(te,1); x_r = results{k}(te,2); 
            plot([x_l+0.1 x_r-0.1]+(k-1)*(n_state+1), y_range(1)+[y_height y_height].*0.57+ast_y_arr(te)*y_height*0.03, 'k-');
            if results{k}(te,6) <= 0.001; ast = ast_arr{3};
            elseif results{k}(te,6) <= 0.01; ast = ast_arr{2};
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

gap_margin = 0.184;
legend_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
for k = 1:K
    ah = axes('position',get(gca,'position'),'visible','off');
    legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
        'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');
    
    ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end

k = K+1;
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah, pl(k), legend_label{k}, 'FontSize',14, ...
    'Location',[0.107+(k-1)*gap_margin 0.905 0.07 0.07], 'Units','normalized');

ax_in = axes(fig,'Position',[0.103+(k-1)*gap_margin 0.68 0.08 0.22], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);

tic;
exportgraphics(fig, sprintf('si_figures/png/mode_dwell_time_w_noto_%s_violin_alpha_20ms_weighted.png',noto_thr), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/mode_dwell_time_w_noto_%s_violin_alpha_20ms_weighted.eps',noto_thr), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/mode_dwell_time_w_noto_%s_violin_alpha_20ms_weighted.fig',noto_thr));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/mode_dwell_time_w_noto_%s_violin_alpha_20ms_weighted.pdf',noto_thr), '-pdf', '-opengl', '-r450');
toc;

%% draw bar plot across bands (probability density)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end
load('../../../results_yjp/UM_info_ch128.mat', 'bands');
explained_all = nan([5,2]);

for fi = 2:7
    band_name = sprintf('band_[%d-%d]',bands(fi,1),bands(fi,2));
    load(sprintf(['%s/tmp_folder/250212_UofM_pca_results/' ...
        '%s/pca_result_of_EC_explained.mat'],dropbox_path,band_name));
    explained_all(fi-1,:) = explained(1:2);
end

fig = figure(1);
fig.Position = [100 100 1200 500]/fs_factor;
clf;

title_fs = 13;
tick_fs = 9;
label_fs = 13;
legend_fs = 12;

% -----------------------------------------------------------------
ax = subplot_tight(2,5,[1,2,6,7],[0.152,0.06]);
hold on; grid on; box on;
x_ticks = 1:6;
x_tick_labels = {'(delta)','(theta)','(alpha)','(beta)_{low}','(beta)_{high}','(gamma)'};
plot(x_ticks,sum(explained_all,2), 'o-', 'Color','#004040', ...
    'LineWidth',3, 'MarkerSize',8, 'MarkerFaceColor','auto');
set(gca, 'FontSize',tick_fs);
xticks(x_ticks);
xticklabels(x_tick_labels);
xlim([0 7]);
ylim([60 90]);
xlabel('Band', 'FontSize',label_fs);
ylabel('Explainability (%)', 'FontSize',label_fs);
legend('(1st PC) + (2nd PC)', 'FontSize',legend_fs);

axes('Position',[.085 .56 .12 .18]);
hold on; grid on; box on;
plot(x_ticks,explained_all(:,1), 'o-', 'Color','#802922', ...
    'LineWidth',2);
xticks(x_ticks);
xticklabels('');
xlim([0.5 6.5]);
ylim([35 60]);
xlabel('Band');
title('1st principal component');

axes('Position',[.24 .56 .12 .18]);
hold on; grid on; box on;
plot(x_ticks,explained_all(:,2), 'o-', 'Color','#BFBE9C', ...
    'LineWidth',2);
xticks(x_ticks);
xticklabels('');
xlim([0.5 6.5]);
ylim([20 27]);
xlabel('Band');
title('2nd principal component');

text(ax, -0.12,1.06, '(A) PCA explainablity over different frequency bands', ...
    'FontSize',title_fs, 'FontWeight','bold', 'Units','normalized');

% -----------------------------------------------------------------
text(ax, 1.06,1.06, '(B) Occurrence for each mode over different frequency bands', ...
    'FontSize',title_fs, 'FontWeight','bold', 'Units','normalized');

K = 4;
state_list = {'EC','LOC','BS','DS','ROC'};
state_label_list = {'Eyes-closed','LOC','Burst','Suppression','ROC'};
n_state = length(state_list);
ls_arr = {'d-','s-','o-','v-','^-','h-'};
occur_all = cell([1,6]);
load_path = sprintf('%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref', ...
    dropbox_path);

mode_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;

for fi = 2:7
    band_name = sprintf('band_[%d-%d]',bands(fi,1),bands(fi,2));
    if fi<4; resol = 100;
    else; resol = 20; end
    load(sprintf('%s/%s/%dms/clustering_w_1min/properties_all.mat', ...
        load_path,band_name,resol));
    d_size = size(occur);
    occur_all{fi-1} = permute(reshape(occur,[d_size(1),d_size(2)*d_size(3),d_size(4)]),[1,3,2]);
end

subplot_idx = [3 8 9 4 5];
for k = 1:K+1
    subplot_tight(2,5,subplot_idx(k),[0.15 0.04]);
    hold on; grid on; box on;
    for fi = 1:6
        if fi==3; lw=3; a=0; else; lw=1; a=fi; end
        y = mean(occur_all{fi}(k,:,:)*100,3,'omitnan');
%         se = std(occur_all{fi}(k,:,:)*100,0,3,'omitnan')./sqrt(sum(~isnan(occur_all{fi}(k,:,:)*100),3));        
%         errorbar(1:5, y, se, 'k', ...
%         'LineStyle','none', 'CapSize',6, 'Color',[c_list(k,:) 1-0.1*a]);
        p(fi) = plot(1:5, y, ...
            ls_arr{fi}, 'Color',[c_list(k,:) 1-0.1*a], 'LineWidth',lw);
    end
    if k~=5; legend_loc = 'southwest';
    else; legend_loc = 'northwest';
    end
    legend(p(3), mode_label{k}, ...
        'FontSize',tick_fs, 'Location',legend_loc);
    xticks(1:n_state);
    xticklabels(state_label_list);
    xlim([0.5 n_state+0.5]);
%     ylim([0 70]);
    if (k==1) || (k== 4); ylim([10 40]);
    elseif (k==2) || (k==3); ylim([0 25]);
    elseif k==5; ylim([0 70]);
    end
    xlabel('State'); ylabel('PDF (%)');
end

% -----------------------------------------------------------------
ax1 = subplot_tight(2,5,10,[0.16 0.04]);
hold on; grid on;
for fi = 1:6
    if fi==3; lw=3; else; lw=1; end
    plot(fi, ls_arr{fi}, 'Color','k', 'LineWidth',lw);
end
axis off;
xlim([-1 0]);
legend(ax1, x_tick_labels, 'NumColumns',2, 'Position',[0.835 0.3 0.1 0.1]);

% -----------------------------------------------------------------
load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', 'topo_idx_comb');
C_comb.EOEC = permute(C_comb.EOEC,[2,3,1]);
load('references/topoplot_parameters.mat');
gap_margin = 0.04;
for k = 1:K
    ax_in = axes(fig,'Position',[0.785+(k-1)*gap_margin 0.08 0.035 0.2], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    text(ax_in,0.5,1.18,mode_label{k}, 'Units','normalized', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
k = K+1;
ax_in = axes(fig,'Position',[0.785+(k-1)*gap_margin 0.08 0.035 0.2], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);
text(ax_in,0.5,1.18,mode_label{k}, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

tic;
exportgraphics(fig, 'si_figures/png/mode_distribution_across_bands.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/mode_distribution_across_bands.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/mode_distribution_across_bands.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/mode_distribution_across_bands.pdf', '-pdf', '-opengl', '-r450');
toc;

%% draw bar plot across bands (dwell time)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end
load('../../../results_yjp/UM_info_ch128.mat', 'bands');
explained_all = nan([5,2]);

for fi = 2:7
    band_name = sprintf('band_[%d-%d]',bands(fi,1),bands(fi,2));
    load(sprintf(['%s/tmp_folder/250212_UofM_pca_results/' ...
        '%s/pca_result_of_EC_explained.mat'],dropbox_path,band_name));
    explained_all(fi-1,:) = explained(1:2);
end

fig = figure(1);
fig.Position = [100 100 1200 500]/fs_factor;
clf;

title_fs = 13;
tick_fs = 9;
label_fs = 13;
legend_fs = 12;

% -----------------------------------------------------------------
ax = subplot_tight(2,5,[1,2,6,7],[0.152,0.06]);
hold on; box on;
x_ticks = 1:6;
x_tick_labels = {'(delta)','(theta)','(alpha)','(beta)_{low}','(beta)_{high}','(gamma)'};
load('../../../results_yjp/UM_info_ch128.mat','bands');
load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', '*idx*');
C_comb_EOEC = C_comb.EOEC(:,topo_idx_comb);
mode_label = {'Mode1','Mode2','Mode3','Mode4','Other'};
K = 4;
corr_mat = nan([K,6]);
for fi = 2:7
    band_name = sprintf('band_[%d-%d]', bands(fi,1),bands(fi,2));
    load(sprintf(['%s/tmp_folder/250212_UofM_pca_results/' ...
        '%s/pca_mask_result_of_EC.mat'],dropbox_path,band_name));
    corr_mat(:,fi-1) = diag(corr(C_comb_EOEC',mask'));
end
for i = 1:3; xline(0.25*i); end
for i = 1:5; yline(1/6*i); end
xlim([0,1]);
ylim([0,1]);
set(gca, 'FontSize',label_fs);
xticks(0.125:0.25:1); xticklabels(mode_label);
yticks(1/12:1/6:1); yticklabels(x_tick_labels);
ax.XAxisLocation = 'top';
% ax.InnerPosition = [0.06 0.152 0.316 0.696];
ax.InnerPosition = [0.1 0.152 0.27 0.6];
ax.YDir = 'reverse';

h_margin = 0.167; w_margin = 0.248;
for fi = 1:6
    for k = 1:K
    text(ax, 0.047+w_margin*(k-1),0.079+h_margin*(fi-1), ...
        sprintf('%.4f',corr_mat(k,fi)), 'FontSize',label_fs, ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');
    end
end


text(ax, -0.28,1.22, '(C) Correlation between centroids across bands', ...
    'FontSize',title_fs, 'FontWeight','bold', 'Units','normalized');

% -----------------------------------------------------------------
text(ax, 1.08,1.22, '(D) Dwell time for each mode over different frequency bands', ...
    'FontSize',title_fs, 'FontWeight','bold', 'Units','normalized');

K = 4;
state_list = {'EC','LOC','BS','DS','ROC'};
state_label_list = {'Eyes-closed','LOC','Burst','Suppression','ROC'};
n_state = length(state_list);
ls_arr = {'d-','s-','o-','v-','^-','h-'};
dwell_all = cell([1,6]);
load_path = sprintf('%s/Moonbrainlab Raw Data/Michigan McDonnell data/preprocessed_wo_reref', ...
    dropbox_path);

c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;

for fi = 2:7
    band_name = sprintf('band_[%d-%d]',bands(fi,1),bands(fi,2));
    if fi<4; resol = 100;
    else; resol = 20; end
    load(sprintf('%s/%s/%dms/clustering_w_1min/properties_all.mat', ...
        load_path,band_name,resol));
    d_size = size(dwell);
    dwell_all{fi-1} = permute(reshape(dwell,[d_size(1),d_size(2)*d_size(3),d_size(4)]),[1,3,2]);
end

subplot_idx = [3 8 9 4 5];
for k = 1:K+1
    ax = subplot_tight(2,5,subplot_idx(k),[0.15 0.04]);
    hold on; grid on; box on;
    for fi = 1:6
        if fi==3; lw=3; a=0; else; lw=1; a=fi; end
        p(fi) = plot(1:5, squeeze(mean(dwell_all{fi}(k,:,:),3,'omitnan')), ...
            ls_arr{fi}, 'Color',[c_list(k,:) 1-0.1*a], 'LineWidth',lw);
    end
    legend(p(3), mode_label{k}, ...
        'FontSize',label_fs, 'Location','northwest');
    xticks(1:n_state);
    xticklabels(state_label_list);
    xlim([0.5 n_state+0.5]);
    ylim([0 0.45]);
    xlabel('State'); ylabel('Dwell time (s)');
end

% -----------------------------------------------------------------
ax = subplot_tight(2,5,10,[0.16 0.04]);
hold on; grid on;
for fi = 1:6
    if fi==3; lw=3; else; lw=1; end
    plot(fi, ls_arr{fi}, 'Color','k', 'LineWidth',lw);
end
axis off;
xlim([-1 0]);
legend(ax, x_tick_labels, 'NumColumns',2, 'Position',[0.835 0.3 0.1 0.1]);

% -----------------------------------------------------------------
load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat', ...
    'C_comb', 'topo_idx_comb');
C_comb.EOEC = permute(C_comb.EOEC,[2,3,1]);
load('references/topoplot_parameters.mat');
gap_margin = 0.04;
for k = 1:K
    ax_in = axes(fig,'Position',[0.785+(k-1)*gap_margin 0.08 0.035 0.2], 'Units','normalized');
    hold on;
    topo = C_comb.EOEC(:,:,k);  
    topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
        'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    text(ax_in,0.5,1.18,mode_label{k}, 'Units','normalized', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
k = K+1;
ax_in = axes(fig,'Position',[0.785+(k-1)*gap_margin 0.08 0.035 0.2], 'Units','normalized');
hold on;
topo = nan(200);
topo(topo_idx_comb) = normrnd(0,eps,[1,length(topo_idx_comb)]);
topoplot_figure(topo, borderCoords, xx, yy, Coord, ...
    'scatter',0,'line',0);
colormap(slanCM('jet'));
clim([-0.75 0.75]);
text(ax_in,0.5,1.18,mode_label{k}, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

tic;
exportgraphics(fig, 'si_figures/png/mode_dwell_time_across_bands.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/mode_dwell_time_across_bands.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/mode_dwell_time_across_bands.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/mode_dwell_time_across_bands.pdf', '-pdf', '-opengl', '-r450');
toc;

%% compare centroids
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

load('../../references/comb_centroids_20240311/combined_centroids_20240311.mat');
C_comb_EOEC = permute(C_comb.EOEC,[2,3,1]);
load('references/topoplot_parameters.mat');
K = 4;

C_HBN_EC = load('../../references/centroid_data/HBN_control_ec.mat','C').C;
C_HBN_EO = load('../../references/centroid_data/HBN_control_eo.mat','C').C;
C_HBN_EOEC = 0.5.*(C_HBN_EC + C_HBN_EO);
C_HBN_EC_20 = load(sprintf('%s/Moonbrainlab Raw Data/Healthy Brain Network data/HBN_pca_across_bands_250211_yj_to_yh/alpha/pca_mask_result_of_EC_tm10_tw10_sm20_re.mat', ...
    dropbox_path),'mask').mask;

C_UofM_EC = load('../../references/centroid_data/combined_centroids_codes/20240223_wo_reref/k_means_EC_all.mat','C').C;
C_UofM_EO = load('../../references/centroid_data/combined_centroids_codes/20240223_wo_reref/k_means_EO_all.mat','C').C;
C_UofM_EOEC = 0.5.*(C_UofM_EC + C_UofM_EO);
C_UofM_EC_20 = load(sprintf('%s/tmp_folder/250212_UofM_pca_results/band_[8-12]/pca_mask_result_of_EC.mat', ...
    dropbox_path),'mask').mask;

fig = figure(1);
fig.Position = [1000 100 1000 500]/fs_factor;
clf; axis off;
ax = subplot(1,1,1);

text_fs = 10;
corr_fs = 8;
% ==================== dataset robustness ====================
% ---------------- Combined EOEC, UofM, HBN -----------------------
x_pos = 0.02; y_pos = 0.71;
w_space = 0.11; h_space = 0.34; h_margin = 0.38;
f_width = 0.2; f_height = 0.22;
% ---------------- Combined centroids (EO & EC) --------------------
text(ax,-0.09,1.05, '(A) Combined centroids (EO & EC)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos f_width f_height]);
    hold on;
    topo = C_comb_EOEC(:,:,k);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
end
% ---------------- UofM centroids (EO & EC) ------------------------
text(ax,-0.09,1.05-h_margin, '(B) Unversity of Michigan (EO & EC)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos-h_space f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_UM) = C_UofM_EOEC(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    rho = corr(C_comb.EOEC(k,topo_idx_comb)',C_UofM_EOEC(k,topo_idx_UM2comb)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end
% ---------------- HBN centroids (EO & EC) ------------------------
text(ax,-0.09,1.05-2*h_margin-0.04, '(C) Healthy Brain Network (EO & EC)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos-2*h_space f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_HBN) = C_HBN_EOEC(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    rho = corr(C_comb.EOEC(k,topo_idx_comb)',C_HBN_EOEC(k,topo_idx_HBN2comb)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end

% ==================== time resolution robustness ====================
x_pos = 0.54; y_pos = 0.74;
w_space = 0.08; h_space = 0.245; h_margin = 0.305;
f_width = 0.06; f_height = 0.18;
% ---------------- UofM EC centroids (100ms) ------------------------
text(ax,0.515,1.05, '(D) UofM eyes-closed centroids (100ms)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_UM) = C_UofM_EC(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    rho = corr(C_UofM_EOEC(k,:)',C_UofM_EC(k,:)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end
% ---------------- UofM EC centroids (20ms) ------------------------
text(ax,0.515,1.05-h_margin, '(E) UofM eyes-closed principal components (20ms)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos-h_space f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_comb) = C_UofM_EC_20(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    c_max = max(abs(C_UofM_EC_20(k,:)));
    clim(c_max*[-1 1]);
    rho = corr(C_UofM_EC(k,topo_idx_UM2comb)',C_UofM_EC_20(k,:)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end
% ---------------- HBN EC centroids (100ms) ------------------------
text(ax,0.515,1.05-2*h_margin, '(F) HBN eyes-closed centroids (100ms)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos-2*h_space f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_HBN) = C_HBN_EC(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    clim([-0.75 0.75]);
    rho = corr(C_HBN_EOEC(k,:)',C_HBN_EC(k,:)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end
% ---------------- HBN EC centroids (20ms) ------------------------
text(ax,0.515,1.05-3*h_margin, '(G) HBN eyes-closed principal components (20ms)', ...
    'FontSize',text_fs, 'FontWeight','bold');
for k = 1:K
    axes('Position',[x_pos+(k-1)*w_space y_pos-3*h_space f_width f_height]);
    hold on;
    topo = nan(200); topo(topo_idx_comb) = C_HBN_EC_20(k,:);
    topoplot_figure(topo,borderCoords,xx,yy,Coord,'scatter',0,'line',0);
    colormap(slanCM('jet'));
    c_max = max(abs(C_HBN_EC_20(k,:)));
    clim(c_max*[-1 1]);
    rho = corr(C_HBN_EC(k,topo_idx_HBN2comb)',C_HBN_EC_20(k,:)');
    title(sprintf('\\rho=%.4f', rho),...
        'FontSize',corr_fs, 'FontWeight','normal');
end

tic;
exportgraphics(fig, 'si_figures/png/compare_centroids.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/compare_centroids.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/compare_centroids.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/compare_centroids.pdf', '-pdf', '-opengl', '-r450');
toc;

%% save mode_distribution with noto
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

K = 4;
K_list = 1:K+1;
state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
state_label_list = {'Eyes Closed', 'LOC', 'Busrt', 'Suppression', 'ROC'};
n_state = length(state_list);

load('../../../results_yjp/UM_info_ch128.mat');

readme_axis = {'K', 'state', 'subject', 'slice', 'measure'};
readme_last_column = {'1: IDX', ...
    '2: avg. weight for each k', ...
    '3: sum weight for each beta axis', ...
    '4: weighted dwelling time', ...
    '5: switches', ...
    '6: occurences', ...
    '7: dwelling time'};
readme_last_column_trans = {'1: trans_prob_mat', ...
    '2: trans_prob_mat_weighted'};
n_measure = length(readme_last_column); 
n_slice = 5; % 5min -> 1min * 5
resol = 100; % 100ms
unit = resol/1e3;
len_slice = 60/unit;

load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/' ...
    'preprocessed_wo_reref/band_[8-12]/%dms/clustering/'], ...
    dropbox_path,resol);

% thr_ratio = 0.15;
for thr_ratio = 0.05:0.05:0.3
    tic;
    mode_dist_all = nan([K+1,n_state,n_sub,n_slice,n_measure]);
    trans_dist_all = nan([K+1,K+1,n_state,n_sub,n_slice,2]);
    
    thr_arr = nan([1,n_sub]); 
    norm_beta_factor = nan([1,n_sub]); 
    for idx = 1:n_sub
        load([load_path '../regression_mae/regr_st_EC_' subjname{idx} '.mat']);
        thr_arr(idx) = quantile(epsilon,1-thr_ratio);
        norm_beta_factor(idx) = mean(sqrt(beta(1,:).^2+beta(2,:).^2));

        for st = 1:n_state
            state = state_list{st};
            try
                load([load_path '../regression_mae/regr_st_' state '_' subjname{idx} '.mat']);                
                load([load_path 'regr_st_' state '_' subjname{idx} '_n_nulls_1000.mat'], 'IDX');
                IDX(epsilon>thr_arr(idx)) = K+1;
                for sl = 1:n_slice
                    idx_list = 1+(sl-1)*len_slice:sl*len_slice;
                    prop = cal_transition_prop_v2(IDX(idx_list), K+1, unit);
                    [~, w_prop] = cal_regression_clustering(beta(:,idx_list),K+1,unit,0,IDX(idx_list));
        
                    mode_dist_all(:,st,idx,sl,1) = histcounts(IDX(idx_list),0.5:K+1.5,'Normalization','pdf');
                    mode_dist_all(:,st,idx,sl,2) = w_prop.avg_weight./norm_beta_factor(idx);
                    mode_dist_all(:,st,idx,sl,3) = w_prop.sum_weight./norm_beta_factor(idx);
                    mode_dist_all(:,st,idx,sl,4) = w_prop.dwell_time./norm_beta_factor(idx);
                    mode_dist_all(:,st,idx,sl,5) = prop.switches;
                    mode_dist_all(:,st,idx,sl,6) = prop.occurrence;
                    mode_dist_all(:,st,idx,sl,7) = prop.dwell_time;

                    trans_dist_all(:,:,st,idx,sl,1) = prop.trans_mat./sum(prop.trans_mat,1);
                    trans_dist_all(:,:,st,idx,sl,2) = w_prop.trans_prob_mat;
                end
            end
        end
    end
    
    save(sprintf('references/mode_dist_profiles_w_noto_%s_%dms_250306.mat', ...
        strrep(num2str(thr_ratio),'.','p'),resol), ...
        'mode_dist_all', 'trans_dist_all', 'K_list', 'state_list', 'state_label_list', ...
        'readme_last_column', 'readme_axis', 'readme_last_column_trans', ...
        'thr_arr', 'thr_ratio', 'norm_beta_factor', '-v7.3');
    toc;

    % to send JYM =======================================================
    desc_columns = {'K', 'n_slice', 'n_sub', 'n_state'};

    mode_dist_all = permute(mode_dist_all,[1,4,3,2,5]);
    occur = mode_dist_all(:,:,:,:,6);
    occur_w = mode_dist_all(:,:,:,:,2);
    dwell = mode_dist_all(:,:,:,:,7);
    dwell_w = mode_dist_all(:,:,:,:,4);
    trans = trans_dist_all(:,:,:,:,:,1);
    trans_w = trans_dist_all(:,:,:,:,:,2);

    save(sprintf('references/properties_all_alpha_w_noto_%s_%dms.mat', ...
        strrep(num2str(thr_ratio),'.','p'),resol), ...
        'desc_columns', 'state_list', 'n_state', ...
        'occur', 'occur_w', 'dwell', 'dwell_w', 'trans', 'trans_w', '-v7.3');
    % ===================================================================
end

%% save mode_distribution without noto
clear, close all;
if ismac
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    dropbox_path = '/Dropbox';
end

K = 4;
K_list = 1:K;
state_list = {'EC', 'LOC', 'BS', 'DS', 'ROC'};
state_label_list = {'Eyes Closed', 'LOC', 'Busrt', 'Suppression', 'ROC'};
n_state = length(state_list);

load('../../../results_yjp/UM_info_ch128.mat');

readme_axis = {'K', 'state', 'subject', 'slice', 'measure'};
readme_last_column = {'1: IDX', ...
    '2: avg. weight for each k', ...
    '3: sum weight for each beta axis', ...
    '4: weighted dwelling time', ...
    '5: switches', ...
    '6: occurences', ...
    '7: dwelling time'};
readme_last_column_trans = {'1: trans_prob_mat', ...
    '2: trans_prob_mat_weighted'};
n_measure = length(readme_last_column); 
n_slice = 5; % 5min -> 1min * 5
resol = 100; % 100ms
unit = resol/1e3;
len_slice = 60/unit;

load_path = sprintf(['%s/Moonbrainlab Raw Data/Michigan McDonnell data/' ...
    'preprocessed_wo_reref/band_[8-12]/%dms/clustering/'], ...
    dropbox_path,resol);

tic;
mode_dist_all = nan([K,n_state,n_sub,n_slice,n_measure]);
trans_dist_all = nan([K,K,n_state,n_sub,n_slice,2]);

thr_arr = nan([1,n_sub]); 
norm_beta_factor = nan([1,n_sub]); 
for idx = 1:n_sub
    for st = 1:n_state
        state = state_list{st};
        try
            load([load_path '../regression_mae/regr_st_' state '_' subjname{idx} '.mat']);                
            load([load_path 'regr_st_' state '_' subjname{idx} '_n_nulls_1000.mat'], 'IDX');
            for sl = 1:n_slice
                idx_list = 1+(sl-1)*len_slice:sl*len_slice;
                prop = cal_transition_prop_v2(IDX(idx_list), K, unit);
                [~, w_prop] = cal_regression_clustering(beta(:,idx_list),K,unit,0,IDX(idx_list));
    
                mode_dist_all(:,st,idx,sl,1) = histcounts(IDX(idx_list),0.5:K+0.5,'Normalization','pdf');
                mode_dist_all(:,st,idx,sl,2) = w_prop.avg_weight./norm_beta_factor(idx);
                mode_dist_all(:,st,idx,sl,3) = w_prop.sum_weight./norm_beta_factor(idx);
                mode_dist_all(:,st,idx,sl,4) = w_prop.dwell_time./norm_beta_factor(idx);
                mode_dist_all(:,st,idx,sl,5) = prop.switches;
                mode_dist_all(:,st,idx,sl,6) = prop.occurrence;
                mode_dist_all(:,st,idx,sl,7) = prop.dwell_time;

                trans_dist_all(:,:,st,idx,sl,1) = prop.trans_mat./sum(prop.trans_mat,1);
                trans_dist_all(:,:,st,idx,sl,2) = w_prop.trans_prob_mat;
            end
        end
    end
end
    
save(sprintf('references/mode_dist_profiles_wo_noto_%dms_250306.mat',resol), ...
    'mode_dist_all', 'trans_dist_all', 'K_list', 'state_list', 'state_label_list', ...
    'readme_last_column', 'readme_axis', 'readme_last_column_trans', ...
    'norm_beta_factor', '-v7.3');
toc;

% to send JYM =======================================================
desc_columns = {'K', 'n_slice', 'n_sub', 'n_state'};

mode_dist_all = permute(mode_dist_all,[1,4,3,2,5]);
occur = mode_dist_all(:,:,:,:,6);
occur_w = mode_dist_all(:,:,:,:,2);
dwell = mode_dist_all(:,:,:,:,7);
dwell_w = mode_dist_all(:,:,:,:,4);
trans = trans_dist_all(:,:,:,:,:,1);
trans_w = trans_dist_all(:,:,:,:,:,2);

save(sprintf('references/properties_all_alpha_wo_noto_%dms.mat',resol), ...
    'desc_columns', 'state_list', 'n_state', ...
    'occur', 'occur_w', 'dwell', 'dwell_w', 'trans', 'trans_w', '-v7.3');
% ===================================================================

%% noto analysis
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

K = 4;
c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
state_label_list = {'Eyes-closed', 'LOC', 'Burst', 'Suppression', 'ROC'};
fs_label = 15;
fs_ticks = 12;
fs_legend = 13;
fs_title = 16;

fig = figure(1);
fig.Position = [100 200 1400 800]/fs_factor;
clf;

marker_list = {'^-', 's-', 'o-', 'd-', 'p-'};
h_space = 0.156;
w_space = 0.06;

% ===================================================
subplot_tight(2,1,1,[h_space w_space]);
hold on; grid on; box on;
cnt = 0;
for noto_thr = 0.05:0.05:0.25
    noto_thr_s = strrep(num2str(noto_thr),'.','p');

    load(sprintf('references/properties_all_alpha_%s.mat',noto_thr_s), ...
        'occur', 'n_state');
    test = squeeze(mean(occur,[2,3],'omitnan'))*100;
    cnt = cnt+1;
    if cnt==3; lw = 4;
    else; lw = 1; end
    for k = 1:K+1
        plot((1:n_state)+(n_state+1)*(k-1), test(k,:), marker_list{cnt}, ...
            'Color',c_list(k,:), 'LineWidth',lw, 'MarkerSize',14);
    end
end
set(gca, 'FontSize',fs_ticks);
xlabel('State', 'FontSize',fs_label);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) ...
     (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list ...
    state_label_list state_label_list]);
ylabel('Probability density (%)', 'FontSize',fs_label);
ylim([0 41]);
for k = 1:K+1
    p(2*(k-1)+1) = plot(1:2, 1e3*(1:2), '-', ...
        'Color',c_list(k,:), 'LineWidth',10, ...
        'MarkerEdgeColor',c_list(k,:));
end
for c = 1:5
    p(2*c) = plot(1:2, 1e3*(1:2), marker_list{c}, ...
        'Color',[0 0 0], 'LineWidth',2, 'MarkerSize',10);
end
legend(p, {'Mode1', '5%', ...
    'Mode2', '10%', ...
    'Mode3', '15%', ...
    'Mode4', '20%', ...
    'Other', '25%'}, ...
    'NumColumns',5, 'FontSize',fs_legend, ...
    'Position',[0.22 0.76 0.38 0.06], 'Units','normalized');
text(-0.05,1.07,'(A) Probability density across thresholds', ...
    'FontSize',fs_title, 'FontWeight','bold', 'Units','normalized');

% ===================================================
subplot_tight(2,1,2,[h_space w_space]);
hold on; grid on; box on;
cnt = 0;
for noto_thr = 0.05:0.05:0.25
    noto_thr_s = strrep(num2str(noto_thr),'.','p');

    load(sprintf('references/properties_all_alpha_%s.mat',noto_thr_s), ...
        'dwell', 'n_state');
    test = squeeze(mean(dwell,[2,3],'omitnan'));
    cnt = cnt+1;
    if cnt==3; lw = 4;
    else; lw = 1; end
    for k = 1:K+1
        plot((1:n_state)+(n_state+1)*(k-1), test(k,:), marker_list{cnt}, ...
            'Color',c_list(k,:), 'LineWidth',lw, 'MarkerSize',14);
    end
end
set(gca, 'FontSize',fs_ticks);
xlabel('State', 'FontSize',fs_label);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) ...
     (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list ...
    state_label_list state_label_list]);
ylabel('Dwell time (s)', 'FontSize',fs_label);
ylim([0.09 0.27]);
for k = 1:K+1
    p(2*(k-1)+1) = plot(1:2, 1e3*(1:2), '-', ...
        'Color',c_list(k,:), 'LineWidth',10, ...
        'MarkerEdgeColor',c_list(k,:));
end
for c = 1:5
    p(2*c) = plot(1:2, 1e3*(1:2), marker_list{c}, ...
        'Color',[0 0 0], 'LineWidth',2, 'MarkerSize',10);
end
legend(p, {'Mode1', '5%', ...
    'Mode2', '10%', ...
    'Mode3', '15%', ...
    'Mode4', '20%', ...
    'Other', '25%'}, ...
    'NumColumns',5, 'FontSize',fs_legend, ...
    'Position',[0.22 0.34 0.38 0.06], 'Units','normalized');
text(-0.05,1.07,'(B) Dwell time across thresholds', ...
    'FontSize',fs_title, 'FontWeight','bold', 'Units','normalized');

tic;
exportgraphics(fig, 'si_figures/png/noto_trends_unweighted.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/noto_trends_unweighted.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/noto_trends_unweighted.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/noto_trends_unweighted.pdf', '-pdf', '-opengl', '-r450');
toc;

%% noto analysis with weigthed
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

K = 4;
c_list = [[0 113 193]; [126 46 141]; [163 20 48]; [237 177 31]; [119 172 48]];
c_list = c_list/256;
state_label_list = {'Eyes-closed', 'LOC', 'Burst', 'Suppression', 'ROC'};
fs_label = 15;
fs_ticks = 12;
fs_legend = 13;
fs_title = 16;

fig = figure(1);
fig.Position = [100 200 1400 800]/fs_factor;
clf;

marker_list = {'^-', 's-', 'o-', 'd-', 'p-'};
h_space = 0.156;
w_space = 0.06;

% ===================================================
subplot_tight(2,1,1,[h_space w_space]);
hold on; grid on; box on;
cnt = 0;
for noto_thr = 0.05:0.05:0.25
    noto_thr_s = strrep(num2str(noto_thr),'.','p');

    load(sprintf('references/properties_all_alpha_%s.mat',noto_thr_s), ...
        'occur_w', 'n_state');
    test = squeeze(mean(occur_w,[2,3],'omitnan'))*100;
    cnt = cnt+1;
    if cnt==2; lw = 4;
    else; lw = 1; end
    for k = 1:K+1
        plot((1:n_state)+(n_state+1)*(k-1), test(k,:), marker_list{cnt}, ...
            'Color',c_list(k,:), 'LineWidth',lw, 'MarkerSize',14);
    end
end
set(gca, 'FontSize',fs_ticks);
xlabel('State', 'FontSize',fs_label);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) ...
     (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list ...
    state_label_list state_label_list]);
ylabel(sprintf('Normalized\noccurrence (%%)'), 'FontSize',fs_label);
ylim([0 48]);
for k = 1:K+1
    p(2*(k-1)+1) = plot(1:2, 1e3*(1:2), '-', ...
        'Color',c_list(k,:), 'LineWidth',10, ...
        'MarkerEdgeColor',c_list(k,:));
end
for c = 1:5
    p(2*c) = plot(1:2, 1e3*(1:2), marker_list{c}, ...
        'Color',[0 0 0], 'LineWidth',2, 'MarkerSize',10);
end
legend(p, {'Mode1', '5%', ...
    'Mode2', '10%', ...
    'Mode3', '15%', ...
    'Mode4', '20%', ...
    'Other', '25%'}, ...
    'NumColumns',5, 'FontSize',fs_legend, ...
    'Position',[0.22 0.76 0.38 0.06], 'Units','normalized');
text(-0.05,1.07,'(A) Normalized occurrence across thresholds', ...
    'FontSize',fs_title, 'FontWeight','bold', 'Units','normalized');

% ===================================================
subplot_tight(2,1,2,[h_space w_space]);
hold on; grid on; box on;
cnt = 0;
for noto_thr = 0.05:0.05:0.25
    noto_thr_s = strrep(num2str(noto_thr),'.','p');

    load(sprintf('references/properties_all_alpha_%s.mat',noto_thr_s), ...
        'dwell_w', 'n_state');
    test = squeeze(mean(dwell_w,[2,3],'omitnan'));
    cnt = cnt+1;
    if cnt==2; lw = 4;
    else; lw = 1; end
    for k = 1:K+1
        plot((1:n_state)+(n_state+1)*(k-1), test(k,:), marker_list{cnt}, ...
            'Color',c_list(k,:), 'LineWidth',lw, 'MarkerSize',14);
    end
end
set(gca, 'FontSize',fs_ticks);
xlabel('State', 'FontSize',fs_label);
xticks([1:n_state (1:n_state)+(n_state+1) (1:n_state)+2*(n_state+1) ...
     (1:n_state)+3*(n_state+1) (1:n_state)+4*(n_state+1)]);
xticklabels([state_label_list state_label_list state_label_list ...
    state_label_list state_label_list]);
ylabel(sprintf('Normalized\ndwell time (s)'), 'FontSize',fs_label);
ylim([0 0.34]);
for k = 1:K+1
    p(2*(k-1)+1) = plot(1:2, 1e3*(1:2), '-', ...
        'Color',c_list(k,:), 'LineWidth',10, ...
        'MarkerEdgeColor',c_list(k,:));
end
for c = 1:5
    p(2*c) = plot(1:2, 1e3*(1:2), marker_list{c}, ...
        'Color',[0 0 0], 'LineWidth',2, 'MarkerSize',10);
end
legend(p, {'Mode1', '5%', ...
    'Mode2', '10%', ...
    'Mode3', '15%', ...
    'Mode4', '20%', ...
    'Other', '25%'}, ...
    'NumColumns',5, 'FontSize',fs_legend, ...
    'Position',[0.22 0.34 0.38 0.06], 'Units','normalized');
text(-0.05,1.07,'(B) Normalized Dwell time across thresholds', ...
    'FontSize',fs_title, 'FontWeight','bold', 'Units','normalized');

tic;
exportgraphics(fig, 'si_figures/png/noto_trends_weighted.png', 'Resolution',450);
exportgraphics(fig, 'si_figures/eps/noto_trends_weighted.eps', 'Resolution',450);
savefig(fig, 'si_figures/fig/noto_trends_weighted.fig');
fig.Color = 'w';
export_fig(fig, 'si_figures/pdf/noto_trends_weighted.pdf', '-pdf', '-opengl', '-r450');
toc;

%% transition matrix 20ms / 100ms alpha band (unweighted)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

fs_title = 11;
fs_ticks = 10;
alphabet = 'A':'Z';
w_space = 0.1846; h_space = 0.43;
w_margin = 0.076; h_margin = 0.16;

file_name = 'w_noto_0p15';
% file_name = 'wo_noto';


fig = figure(1);
fig.Position = [100 100 1200 500]/fs_factor;
clf;

load(['references/mode_dist_profiles_' file_name '_100ms_250306.mat']);
test = mean(trans_dist_all,[4,5],'omitnan');
for st = 1:5
    subplot_tight(2,5,st,[h_margin w_margin]);
    h = heatmap(test(:,:,st,:,:,1)');
    h.CellLabelFormat = '%.2f';
    h.Colormap = slanCM('hot');
    h.XLabel = 'Mode (to)';
    h.YLabel = 'Mode (from)';
    h.ColorLimits = [0 1];
    annotation('textbox', [0.038+(st-1)*w_space, 0.8, 0.1, 0.1], ...
        'String',sprintf('(%s) %s (100ms)',alphabet(st),state_label_list{st}), ...
        'EdgeColor','none', 'FontSize',fs_title, 'FontWeight','bold', ...
        'HorizontalAlignment','left');
    set(gca, 'FontSize',fs_ticks);
end

load(['references/mode_dist_profiles_' file_name '_20ms_250306.mat']);
test = mean(trans_dist_all,[4,5],'omitnan');
for st = 1:5
    subplot_tight(2,5,st+5,[h_margin w_margin]);
    h = heatmap(test(:,:,st,:,:,1)');
    h.CellLabelFormat = '%.2f';
    h.Colormap = slanCM('hot');
    h.XLabel = 'Mode (to)';
    h.YLabel = 'Mode (from)';
    h.ColorLimits = [0 1];
    annotation('textbox', [0.038+(st-1)*w_space, 0.8-h_space, 0.1, 0.1], ...
        'String',sprintf('(%s) %s (20ms)',alphabet(st+5),state_label_list{st}), ...
        'EdgeColor','none', 'FontSize',fs_title, 'FontWeight','bold', ...
        'HorizontalAlignment','left');
    set(gca, 'FontSize',fs_ticks);
end

tic;
exportgraphics(fig, sprintf('si_figures/png/transition_matrix_%s_unweighted.png',file_name), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/transition_matrix_%s_unweighted.eps',file_name), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/transition_matrix_%s_unweighted.fig',file_name));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/transition_matrix_%s_unweighted.pdf',file_name), '-pdf', '-opengl', '-r450');
toc;

%% transition matrix 20ms / 100ms alpha band (weighted)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
    dropbox_path = '/Volumes/YJ_extHDD/Dropbox';
elseif ispc
    fs_factor = 1;
    dropbox_path = '/Dropbox';
end

fs_title = 11;
fs_ticks = 10;
alphabet = 'A':'Z';
w_space = 0.1846; h_space = 0.43;
w_margin = 0.076; h_margin = 0.16;

file_name = 'w_noto_0p15';
% file_name = 'wo_noto';


fig = figure(1);
fig.Position = [100 100 1200 500]/fs_factor;
clf;

load(['references/mode_dist_profiles_' file_name '_100ms_250306.mat']);
test = mean(trans_dist_all,[4,5],'omitnan');
for st = 1:5
    subplot_tight(2,5,st,[h_margin w_margin]);
    h = heatmap(test(:,:,st,:,:,2)');
    h.CellLabelFormat = '%.2f';
    h.Colormap = slanCM('hot');
    h.XLabel = 'Mode (to)';
    h.YLabel = 'Mode (from)';
    h.ColorLimits = [0 1];
    annotation('textbox', [0.038+(st-1)*w_space, 0.8, 0.1, 0.1], ...
        'String',sprintf('(%s) %s (100ms)',alphabet(st),state_label_list{st}), ...
        'EdgeColor','none', 'FontSize',fs_title, 'FontWeight','bold', ...
        'HorizontalAlignment','left');
    set(gca, 'FontSize',fs_ticks);
end

load(['references/mode_dist_profiles_' file_name '_20ms_250306.mat']);
test = mean(trans_dist_all,[4,5],'omitnan');
for st = 1:5
    subplot_tight(2,5,st+5,[h_margin w_margin]);
    h = heatmap(test(:,:,st,:,:,2)');
    h.CellLabelFormat = '%.2f';
    h.Colormap = slanCM('hot');
    h.XLabel = 'Mode (to)';
    h.YLabel = 'Mode (from)';
    h.ColorLimits = [0 1];
    annotation('textbox', [0.038+(st-1)*w_space, 0.8-h_space, 0.1, 0.1], ...
        'String',sprintf('(%s) %s (20ms)',alphabet(st+5),state_label_list{st}), ...
        'EdgeColor','none', 'FontSize',fs_title, 'FontWeight','bold', ...
        'HorizontalAlignment','left');
    set(gca, 'FontSize',fs_ticks);
end

tic;
exportgraphics(fig, sprintf('si_figures/png/transition_matrix_%s_weighted.png',file_name), 'Resolution',450);
exportgraphics(fig, sprintf('si_figures/eps/transition_matrix_%s_weighted.eps',file_name), 'Resolution',450);
savefig(fig, sprintf('si_figures/fig/transition_matrix_%s_weighted.fig',file_name));
fig.Color = 'w';
export_fig(fig, sprintf('si_figures/pdf/transition_matrix_%s_weighted.pdf',file_name), '-pdf', '-opengl', '-r450');
toc;

%% tmp

% means = squeeze(mean(occur,[2,3],'omitnan'));
% stds = squeeze(std(occur,0,[2,3],'omitnan')./sqrt(sum(~isnan(occur),[2,3])));

% means = squeeze(mean(dwell,[2,3],'omitnan'));
% stds = squeeze(std(dwell,0,[2,3],'omitnan')./sqrt(sum(~isnan(dwell),[2,3])));

means = squeeze(mean(data,[1,2],'omitnan'));
stds = squeeze(std(data,0,[1,2],'omitnan')./sqrt(sum(~isnan(data),[1,2])));
