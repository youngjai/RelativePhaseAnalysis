%% save all V and beta
clear, close all;

n_sub = 22;
n_Yeo = 17;
n_TR = 97; % except first TR to synchronize EEG timing

load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri';
file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';
Yeo = niftiread(sprintf(['%s/fmri_references/' ...
    'Yeo2011_%dNetworks_MNI152_FreeSurferConformed3mm_LiberalMask.nii'], ...
    load_path,n_Yeo));

V_all = nan([size(Yeo), n_TR, n_sub]);
Yeo_all = nan([n_Yeo, n_TR, n_sub]);
beta_all = nan([n_TR*21, n_sub]);

for sub = 1:n_sub
    disp(sub);
    try
        % BOLD of all voxels
        V = niftiread(sprintf(['%s/preproc/sub-%02d/ses-01/func/' ...
            'sub-%02d_ses-01_task-checker_bold/func_preproc/%s'], ...
            load_path,sub,sub,file_name));
        V_all(:,:,:,:,sub) = V(:,:,:,2:end);

        % BOLD of all Yeo's networks
        V = permute(V,[4,1,2,3]);
        for y = 1:n_Yeo
            Yeo_all(y,:,sub) = mean(V(2:end,Yeo==y),2,'omitnan');
        end

        % EEG beta1
        load(sprintf('../regression/CB_EC/regr_st_CB_EC_subj_%02d.mat',sub), 'beta');
        beta_all(:,sub) = beta(1,:);
    end
end

Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicA', 'N10: LimbicB', 'N11: ContC', 'N12: ContA', ...
    'N13: ContB', 'N14: TempPar', 'N15: DefaultC', 'N16: DefaultA', ...
    'N17: DefaultB'};

Yeo17_label_abbr = {'VC', 'VP', 'SMA', 'SMB', ...
    'DAA', 'DAB', 'VAA', 'VAB', ...
    'LA', 'LB', 'FPC', 'FPA', ...
    'FPB', 'TP', 'DMC', 'DMA', 'DMB'};

save('all_data_in_NKI_w_gsr_checkeron.mat', ...
    'n_sub', 'n_TR', 'n_Yeo', 'Yeo', ...
    'Yeo_all', 'V_all', 'beta_all', ...
    'Yeo17_label', 'Yeo17_label_abbr', '-v7.3');

%% edit schaefer_fsa5 & save
clear, close all;

n_parcel = 1000;
[~, Ll, ctl] = read_annotation(sprintf('../fmri_references/lh.Schaefer2018_%dParcels_17Networks_order.annot',n_parcel));
[~, Lr, ctr] = read_annotation(sprintf('../fmri_references/rh.Schaefer2018_%dParcels_17Networks_order.annot',n_parcel));
for i = 1:n_parcel/2+1
    Ll(Ll==ctl.table(i,5)) = i-1;
end
Lr(Lr==ctr.table(1,5)) = 0;
for i = 2:n_parcel/2+1
    Lr(Lr==ctr.table(i,5)) = i-1+n_parcel/2;
end
label_vector = cat(1,Ll,Lr);
csvwrite(sprintf('schaefer_%d_fsa5_yj.csv',n_parcel), label_vector);


%% mean topo vector of NKI checkeron EEG data (save)
clear, close all;

load('../../mfiles/references/comb_centroids_20240311/combined_centroids_20240311.mat', 'topo_idx_comb');

scale = 21; % duration between two frames : 2.1s/21 = 0.1s
idx_list_EC_eeg = [ 2*scale: 9*scale 21*scale:28*scale 40*scale:47*scale ...
    59*scale:66*scale 78*scale:85*scale]+1*scale;
idx_list_CB_eeg = [11*scale:18*scale 30*scale:37*scale 49*scale:56*scale ...
    68*scale:75*scale 87*scale:94*scale]+1*scale;

topo_vector_all = nan([2037,31739,22]);
for sub = 1:22
    try
        tmp = load(sprintf('../data/topo_20ms/topo_100ms_wo_reref_cb/sub-%02d ses-01 task-checker eeg GAPArm 100ms_topo_wo_reref.mat',sub),'topo').topo;
        T = size(tmp,1);
        for t = 1:T
            topo_vector_all(t,:,sub) = tmp{t}(topo_idx_comb);
        end
    end
end
n_sub_wo_nan = sum(~isnan(topo_vector_all(1,1,:)));
topo_vector_all = mean(topo_vector_all,3,'omitnan');
save('nki_mean_topo_vector_of_checkeron.mat', ...
    'topo_vector_all', 'n_sub_wo_nan', 'idx_list_EC_eeg', 'idx_list_CB_eeg', '-v7.3');


%% draw combined centroid (UofM and HBN) (save)
clear, close all;

load('../../preprocessing/completed_240715/UM_8_wo_badchan_and_wo_reref.mat');

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
save('topoplot_parameters.mat', ...
    'xx', 'yy', 'Coord', 'borderCoords', '-v7.3');

%% draw difference BOLD and EEG images
clear, close all;
% user_addpath(false, false);
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
elseif ispc
    fs_factor = 1;
end

% fMRI image ==============================================================
n_parcel = 1000;
delay = 6;
n_Yeo = 17;

% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
schaefer = niftiread(['references/NKI/fmri_references/Schaefer2018_' num2str(n_parcel) ...
    'Parcels_' num2str(n_Yeo) 'Networks_order_FSLMNI152_2mm.nii.gz']);
n_ROI = length(unique(schaefer));
[nx_tmp,ny_tmp,nz_tmp] = size(schaefer);

load('references/NKI/all_data_in_NKI_w_gsr_checkeron.mat', 'V_all');
V = mean(V_all,5,'omitnan'); V = zscore(V,0,4);
[nx,ny,nz,nT] = size(V);
V_upsample = nan([nT,nx_tmp,ny_tmp,nz_tmp]);
for t = 1:nT
    V_upsample(t,:,:,:) = imresize3(V(:,:,:,t),[nx_tmp,ny_tmp,nz_tmp],'Method','nearest');
end

idx_list_EC = [ 2: 9 21:28 40:47 59:66 78:85]+delay;
idx_list_CB = [11:18 30:37 49:56 68:75 87:94]+delay;

avg_array = zeros([n_ROI,nT]);
for roi = 2:n_ROI
    idx_roi = find(schaefer==roi-1);
    avg_array(roi,:) = mean(V_upsample(:,idx_roi),2,'omitnan');
end
avg_array_z = zscore(avg_array,0,2);

avg_array_fsa5.EC = parcel_to_surface(mean(avg_array_z(:,idx_list_EC(idx_list_EC<=nT)),2), sprintf('references/NKI/schaefer_%d_fsa5_yj',n_parcel));
avg_array_fsa5.CB = parcel_to_surface(mean(avg_array_z(:,idx_list_CB(idx_list_CB<=nT)),2), sprintf('references/NKI/schaefer_%d_fsa5_yj',n_parcel));

diff_fsa5 = avg_array_fsa5.EC - avg_array_fsa5.CB;

% EEG image ===============================================================
delay_eeg = 21*delay;

load('references/NKI/nki_mean_topo_vector_of_checkeron.mat');
load('references/comb_centroids_20240311/combined_centroids_20240311.mat', 'topo_idx_comb');
load('references/topoplot_parameters.mat');

topo_vector_EC = (mean(topo_vector_all(idx_list_EC_eeg(idx_list_EC_eeg<=2037-delay_eeg)+delay_eeg,:),1)-mean(topo_vector_all,1))./std(topo_vector_all,0,1);
topo_vector_CB = (mean(topo_vector_all(idx_list_CB_eeg(idx_list_CB_eeg<=2037-delay_eeg)+delay_eeg,:),1)-mean(topo_vector_all,1))./std(topo_vector_all,0,1);

diff_topo = topo_vector_EC - topo_vector_CB;

% draw image ==============================================================
fig = figure(1);
fig.Position = [100 100 1200 600]/fs_factor;
clf;
pp = plot_cortical(diff_fsa5, 'surface_name','fsa5');

h = 0.22;
w = 0.16;
for i = 1:4
    pp(i).CLim = [-1.0 1.0]; 
    pp(i).Position([1 3 4]) = [0.1+(i-1)*w w h];
end
set(gca, 'FontSize',10);
cb = get(fig, 'Children'); cb = cb(1);
cb.Limits = [-1 1];
cb.Label.String = 'BOLD difference between EC and CB';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'normal';
cb.Location = 'north';
cb.Position = [0.28 0.57 0.28 0.02];
cb.Ticks = [-1 0 1];
cb.TickLabels = {'(CB related)\newline        -1' 0 '(EC related)\newline         1'};

axes('Position',[0.12+i*w 0.31 w 0.28]);
topo = nan(200);
topo(topo_idx_comb) = diff_topo;
topoplot_figure(topo,borderCoords,xx,yy,Coord);
text(0.51,1.15,{'EEG relative phase' 'difference between EC and CB'}, ...
    'Units','normalized', 'FontSize',11, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

colormap(flip(slanCM('PiYG')));
cb = colorbar;
cb.Limits = [-0.25 0.25];
cb.Position = [0.91 0.32 0.01 0.24];
cb.Ticks = [-0.25 0 0.25];
cb.TickLabels = {'     -0.25\newline(CB related)' '         0' ...
    '      0.25\newline(EC related)'};

tic;
exportgraphics(fig, ['main_figures_NKI/png/difference_fMRI_EEG_map_' ...
    num2str(n_parcel) 'Parcels_w_delay_' num2str(delay) 's.png'], 'Resolution',450);
exportgraphics(fig, ['main_figures_NKI/eps/difference_fMRI_EEG_map_' ...
    num2str(n_parcel) 'Parcels_w_delay_' num2str(delay) 's.eps'], 'Resolution',450);
savefig(fig, ['main_figures_NKI/fig/difference_fMRI_EEG_map_' ...
    num2str(n_parcel) 'Parcels_w_delay_' num2str(delay) 's.fig']);
fig.Color = 'w';
export_fig(fig, ['main_figures_NKI/pdf/difference_fMRI_EEG_map_' ...
    num2str(n_parcel) 'Parcels_w_delay_' num2str(delay) 's.pdf'], '-pdf', '-opengl', '-r450');
toc;

%% draw mean trends beta and yeo network
clear, close all;
% user_addpath(false, false);
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
elseif ispc
    fs_factor = 1;
end

load('references/NKI/all_data_in_NKI_w_gsr_checkeron.mat', 'beta_all', 'Yeo_all');

beta_all_TR = moving_time_window_v2(beta_all,21,105,@(x) x);
beta1_z_mean = mean(zscore(beta_all_TR,0,1),2,'omitnan');
Yeo_z_mean = mean(zscore(Yeo_all,0,2),3,'omitnan');

fig = figure(1);
fig.Position = [100 100 1200 700]/fs_factor;
clf;
subplot(3,2,[3 5]);
hold on; grid on;
set(gca, 'FontSize',16);
for bl = 0:4
    b(1) = patch([ 0 20 20  0]+40*bl+6, [-2 -2 2 2], 'w', ...
        'FaceAlpha',0.2, 'LineStyle',':', 'LineWidth',0.5);
    b(2) = patch([20 40 40 20]+40*bl+6, [-2 -2 2 2], 'g', ...
        'FaceAlpha',0.2, 'LineStyle',':', 'LineWidth',0.5);
end
p(1) = plot((1:97)*2.1, Yeo_z_mean(16,:), 'k-', 'LineWidth',5);
p(2) = plot((1:93)*2.1, beta1_z_mean, 'r-', 'LineWidth',5);
ylim([-1.2 1.2]);
xlim([0 200]);
xticks(0:40:200);
xlabel('Time (s)', 'FontSize',18);
ylabel('Activity', 'FontSize',18);
legend(b, {'Eyes-closed resting', 'Checkerboard'}, ...
    'NumColumns',2, 'Location','north', 'FontSize',12);
title(sprintf('Default Mode Network\nvs. Front-to-Back Mode'), ...
    'FontSize',16, 'Position',[100 1.8 0], 'FontWeight','normal');

ax = axes(fig, 'Position',[0.16 0.64 0.27 0.1]);
box on; hold on;
plot([0.26-0.1 0.26+0.1], [0.64 0.64], 'k-', 'LineWidth',4);
plot([1-0.26-0.1 1-0.26+0.1], [0.64 0.64], 'r-', 'LineWidth',4);
xticks([]);
yticks([]);
xlim([0 1]); ylim([0 1]);
text(0.26,0.34,'DMA', 'Units','normalized', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(1-0.26,0.34,'FB mode', 'Units','normalized', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center');

subplot(3,2,[4 6]);
hold on; grid on;
set(gca, 'FontSize',16);
for bl = 0:4
    b(1) = patch([ 0 20 20  0]+40*bl+6, [-2 -2 2 2], 'w', ...
        'FaceAlpha',0.2, 'LineStyle',':', 'LineWidth',0.5);
    b(2) = patch([20 40 40 20]+40*bl+6, [-2 -2 2 2], 'g', ...
        'FaceAlpha',0.2, 'LineStyle',':', 'LineWidth',0.5);
end
plot((1:97)*2.1, Yeo_z_mean(1,:), 'k-.', 'LineWidth',5);
plot((1:93)*2.1, -beta1_z_mean, 'r-.', 'LineWidth',5);
ylim([-1.2 1.2]);
xlim([0 200]);
xticks(0:40:200);
xlabel('Time (s)', 'FontSize',18);
ylabel('Activity', 'FontSize',18);
legend(b, {'Eyes-closed resting', 'Checkerboard'}, ...
    'NumColumns',2, 'Location','north', 'FontSize',12);
title(sprintf('Visual Central Network\nvs. Back-to-Front Mode'), ...
    'FontSize',16, 'Position',[100 1.8 0], 'FontWeight','normal');

ax = axes(fig, 'Position',[0.598 0.64 0.27 0.1]);
box on; hold on;
plot([0.26-0.1 0.26+0.1], [0.64 0.64], 'k-.', 'LineWidth',4.5);
plot([1-0.26-0.1 1-0.26+0.1], [0.64 0.64], 'r-.', 'LineWidth',4.5);
xticks([]);
yticks([]);
xlim([0 1]); ylim([0 1]);
text(0.26,0.34,'VC', 'Units','normalized', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(1-0.26,0.34,'BF mode', 'Units','normalized', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center');

tic;
exportgraphics(fig, 'main_figures_NKI/png/mean_trend_of_beta_BOLD.png', 'Resolution',450);
exportgraphics(fig, 'main_figures_NKI/eps/mean_trend_of_beta_BOLD.eps', 'Resolution',450);
savefig(fig, 'main_figures_NKI/fig/mean_trend_of_beta_BOLD.fig');
fig.Color = 'w';
export_fig(fig, 'main_figures_NKI/pdf/mean_trend_of_beta_BOLD.pdf', '-pdf', '-opengl', '-r450');
toc;

%% correlation bar graph (CB)
clear, close all;
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
elseif ispc
    fs_factor = 1;
end

load('references/NKI/all_data_in_NKI_w_gsr_checkeron.mat', ...
    'beta_all', 'Yeo_all', 'Yeo17_label_abbr', 'n_Yeo');

corr_arr = nan([1,n_Yeo]);
pval_arr = nan([1,n_Yeo]);

beta_all_TR = moving_time_window_v2(beta_all,21,105,@(x) x);
for yeo = 1:n_Yeo
    [corr_arr(yeo), pval_arr(yeo)] = ...
        corr(mean(beta_all_TR,2,'omitnan'), ...
            mean(Yeo_all(yeo,1:93,:),3,'omitnan')');
end
pval_arr = -log10(pval_arr);

% -------------------------- draw Yeo's 17 networks ---------------------
n_parcel = 1000;
n_Yeo = 17;

% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
schaefer = niftiread(['references/NKI/fmri_references/Schaefer2018_' num2str(n_parcel) ...
    'Parcels_' num2str(n_Yeo) 'Networks_order_FSLMNI152_2mm.nii.gz']);
n_ROI = length(unique(schaefer));
[nx_tmp,ny_tmp,nz_tmp] = size(schaefer);

Yeo = niftiread(['references/NKI/fmri_references/Yeo2011_' num2str(n_Yeo) ...
    'Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii']);
Yeo = imresize3(Yeo,[nx_tmp,ny_tmp,nz_tmp],'Method','nearest');

avg_array = zeros([n_ROI,1]);
for roi = 2:n_ROI
    idx_roi = find(schaefer==roi-1);
    avg_array(roi) = mode(Yeo(idx_roi));
end

avg_array_fsa5 = parcel_to_surface(avg_array, sprintf('references/NKI/schaefer_%d_fsa5_yj',n_parcel));

% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering
colormap_yeo = load('references/NKI/fmri_references/17NetworksColors.mat').colors./255;

fig = figure(1);
fig.Position = [100 100 800 750]/fs_factor;
clf;
pp = plot_cortical(avg_array_fsa5, 'surface_name','fsa5');
for ax = 1:4
    pp(ax).Position([2,4]) = [0.8 0.18];
    pp(ax).CLim = [0, n_Yeo];
    pp(ax).Colormap = colormap_yeo;
end
colorbar off;

axes(fig, 'Position',[0.13 0.06 0.775 0.64]);
hold on; grid on; box on;
bar(1:n_Yeo, corr_arr, ...
    'CData',pval_arr, 'FaceColor','flat');
set(gca, 'FontSize',16);
colormap(gca, flip(slanCM('bone')));
ylim([-0.701 0.701]);
% clim([1 3]);
caxis([1 3]);
cb = colorbar;
cb.Label.String = 'p-value';
cb.Label.FontSize = 20;
cb.Ticks = 1:3;
cb.TickLabels = 10.^-(1:3);
xticks(1:n_Yeo);
xticklabels(Yeo17_label_abbr);
xtickangle(90);
xlabel("Yeo's 17 network", 'FontSize',20);
ylabel("Pearson's correlation coefficient", 'FontSize',20);
set(gca, 'XAxisLocation','top');

axes(fig, 'Position',[0.13 0.664 0.7017 0.024],'Units','normalized');
axis off;
for y = 1:n_Yeo
    patch([-0.3 -0.3 0.3 0.3]+y, [0 1 1 0], 'k', ...
        'FaceColor',colormap_yeo(y+1,:));
end
xlim([-0.2 n_Yeo+1.2]);
xticks(1:n_Yeo);

tic;
exportgraphics(fig, 'main_figures_NKI/png/yeo17_correlation.png', 'Resolution',450);
exportgraphics(fig, 'main_figures_NKI/eps/yeo17_correlation.eps', 'Resolution',450);
savefig(fig, 'main_figures_NKI/fig/yeo17_correlation.fig');
fig.Color = 'w';
export_fig(fig, 'main_figures_NKI/pdf/yeo17_correlation.pdf', '-pdf', '-opengl', '-r450');
toc;

%% correlation map
clear, close all;
% user_addpath(false, false);
if ismac
    fs_factor = 96/get(groot, 'ScreenPixelsPerInch');
elseif ispc
    fs_factor = 1;
end

n_parcel = 1000;
n_Yeo = 17;

% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
schaefer = niftiread(['references/NKI/fmri_references/Schaefer2018_' num2str(n_parcel) ...
    'Parcels_' num2str(n_Yeo) 'Networks_order_FSLMNI152_2mm.nii.gz']);
n_ROI = length(unique(schaefer));
[nx_tmp,ny_tmp,nz_tmp] = size(schaefer);

load('references/NKI/all_data_in_NKI_w_gsr_checkeron.mat', 'V_all', 'beta_all');
V = mean(V_all,5,'omitnan'); V = zscore(V,0,4);
[nx,ny,nz,nT] = size(V);
V_upsample = nan([nT,nx_tmp,ny_tmp,nz_tmp]);
for t = 1:nT
    V_upsample(t,:,:,:) = imresize3(V(:,:,:,t),[nx_tmp,ny_tmp,nz_tmp],'Method','nearest');
end

avg_array = zeros([n_ROI,nT]);
for roi = 2:n_ROI
    idx_roi = find(schaefer==roi-1);
    avg_array(roi,:) = mean(V_upsample(:,idx_roi),2,'omitnan');
end
avg_array_z = zscore(avg_array,0,2);

beta_all_TR = moving_time_window_v2(beta_all,21,105,@(x) x);
beta1_mean = mean(beta_all_TR,2,'omitnan');

avg_array_z_beta_corr = zeros([size(avg_array_z,1),1]);
for roi = 2:size(avg_array_z,1)
    avg_array_z_beta_corr(roi) = corr(avg_array_z(roi,1:end-4)',beta1_mean);
end

avg_array_fsa5 = parcel_to_surface(avg_array_z_beta_corr, sprintf('references/NKI/schaefer_%d_fsa5_yj',n_parcel));

% draw image ==============================================================
fig = figure(1);
fig.Position = [100 100 1200 600];
clf;
pp = plot_cortical(avg_array_fsa5, 'surface_name','fsa5');

h = 0.22;
w = 0.16;
for i = 1:4
    pp(i).Position([1 3 4]) = [0.1+(i-1)*w w h];
end
set(gca, 'FontSize',10);
cb = get(fig, 'Children'); cb = cb(1);
cb.Limits = [-0.6 0.6];
cb.Label.String = 'Correlation between BOLD and FB mode';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'normal';
cb.Location = 'north';
cb.Position = [0.28 0.57 0.28 0.02];
cb.Ticks = [-0.6 0 0.6];
cb.TickLabels = {'(BF correlated)\newline          -0.6' ...
    0 '(FB correlated)\newline           0.6'};

load('references/comb_centroids_20240311/combined_centroids_20240311.mat', 'C_comb');
load('references/topoplot_parameters.mat');

axes('Position',[0.12+i*w 0.31 w 0.28]);
topoplot_figure(squeeze(C_comb.EOEC(1,:,:)),borderCoords,xx,yy,Coord);
text(0.51,1.15,{'EEG relative phase' 'FB mode reference'}, ...
    'Units','normalized', 'FontSize',11, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
cb = colorbar;
cb.Limits = [-0.75 0.75];
cb.Position = [0.91 0.32 0.01 0.24];
cb.Ticks = [-0.75 0 0.75];
cb.TickLabels = {' -0.75\newline (Lag)' '    0' ...
    ' 0.75\newline(Lead)'};

for i = 1:4
    pp(i).Colormap = flip(slanCM('RdBu'));
    pp(i).CLim = [-0.6 0.6];
end

tic;
exportgraphics(fig, 'main_figures_NKI/png/bold_beta1_correlation_v2.png', 'Resolution',450);
exportgraphics(fig, 'main_figures_NKI/eps/bold_beta1_correlation_v2.eps', 'Resolution',450);
savefig(fig, 'main_figures_NKI/fig/bold_beta1_correlation_v2.fig');
fig.Color = 'w';
export_fig(fig, 'main_figures_NKI/pdf/bold_beta1_correlation_v2.pdf', '-pdf', '-opengl', '-r450');
toc;
