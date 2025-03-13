%% correlation image map individually
clear, close all;
user_addpath(false, false);
addpath('/moonbrainlab/YJPark/matlab/toolbox/ENIGMA-master/matlab/shared/surfaces');

n_parcel = 1000;
delay = 4;
load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';
file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';

% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
schaefer = niftiread(sprintf('../fmri_references/Schaefer2018_%dParcels_17Networks_order_FSLMNI152_2mm.nii.gz',n_parcel));
n_ROI = length(unique(schaefer));
[nx_tmp,ny_tmp,nz_tmp] = size(schaefer);

load('all_yeo_and_beta_CB.mat');

avg_array_z_beta_corr = zeros([n_ROI,n_sub]);

for sub = 1:n_sub
    try
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
            load_path, sub, sub, file_name));
        [nx,ny,nz,nT] = size(V);
        V_upsample = nan([nT,nx_tmp,ny_tmp,nz_tmp]);
        for t = 1:nT
            V_upsample(t,:,:,:) = imresize3(V(:,:,:,t),[nx_tmp,ny_tmp,nz_tmp],'Method','nearest');
        end
        V_upsample = zscore(V_upsample,0,1);
        
        avg_array = zeros([n_ROI,nT]);
        for roi = 2:n_ROI
            idx_roi = find(schaefer==roi-1);
            avg_array(roi,:) = mean(V_upsample(:,idx_roi),2,'omitnan');
        end
        avg_array_z = zscore(avg_array,0,2);
        
        for roi = 2:n_ROI
            avg_array_z_beta_corr(roi,sub) = corr(avg_array_z(roi,2:94)',beta_all_TR(1,:,sub)');
        end
    end
end
save(sprintf('gsr_BOLD_beta_corr_map_parcel%d.mat',n_parcel), ...
    'avg_array_z_beta_corr', '-v7.3');


%% correlation map
clear, close all;
user_addpath(false, false);
addpath('/moonbrainlab/YJPark/matlab/toolbox/ENIGMA-master/matlab/shared/surfaces');
n_parcel = 1000;
n_sub = 22;

load(sprintf('gsr_BOLD_beta_corr_map_parcel%d.mat',n_parcel));

l_margin = 0.05;
b_margin = 0.15;
h_fig = 0.45;
w_fig = 0.2;
gap_margin_ud = 0.35;
gap_margin_rl = 0.15;
gap_middle = 0.4;
c_lim = 1.8;

fig = figure(1);
fig.Position = [100 100 1000 500];

for sub = 1:n_sub
    try
        avg_array_fsa5 = parcel_to_surface(avg_array_z_beta_corr(:,sub), sprintf('schaefer_%d_fsa5_yj',n_parcel));
        
        clf;
        pp_ec = plot_cortical(avg_array_fsa5, 'surface_name','fsa5');
        for p = pp_ec; caxis(p, [-0.25 0.25]); end
        annotation('textbox',[0.45 0.6 0.1 0.1], 'String',sprintf('sub-%02d',sub), ...
            'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
            'FontSize',20, 'EdgeColor','none');
        annotation('textbox',[0.45 0.06 0.1 0.1], 'String',"Pearson's correlation coefficient", ...
            'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
            'FontSize',14, 'EdgeColor','none');
        
        exportgraphics(fig, sprintf('../figures/checkerboard/gsr_BOLD_beta_corr_map_sub-%02d.png',sub));
    end
end

avg_array_fsa5 = parcel_to_surface(mean(avg_array_z_beta_corr,2,'omitnan'), sprintf('schaefer_%d_fsa5_yj',n_parcel));

clf;
pp_ec = plot_cortical(avg_array_fsa5, 'surface_name','fsa5');
for p = pp_ec; caxis(p, [-0.1 0.1]); end
annotation('textbox',[0.45 0.6 0.1 0.1], 'String','sub-all', ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'FontSize',20, 'EdgeColor','none');
annotation('textbox',[0.45 0.06 0.1 0.1], 'String',"Pearson's correlation coefficient", ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'FontSize',14, 'EdgeColor','none');

exportgraphics(fig, '../figures/checkerboard/gsr_BOLD_beta_corr_map_sub-all.png');

%% mean trends of V1, DMN_B, beta_1
clear, close all;
user_addpath(false, false);
addpath('/moonbrainlab/YJPark/matlab/toolbox/ENIGMA-master/matlab/shared/surfaces');

Yeo = niftiread('../fmri_references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};
load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';
file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';

yeo_all = nan([17,98,22]);
n_sub = 22;
for sub = 1:n_sub
    try
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
            load_path, sub, sub, file_name));
        V = zscore(V,0,4);
        V = permute(V, [4,1,2,3]);
        for yeo = 1:17
            idx = find(Yeo==yeo);
            yeo_all(yeo,:,sub) = mean(V(:,idx),2,'omitnan');
        end
    end
end
yeo_mean = mean(yeo_all,3,'omitnan');

beta1_all = nan([2037,22]);
for sub = 1:n_sub
    try
        load(sprintf('../regression/CB_EC/regr_st_CB_EC_subj_%02d.mat',sub), 'beta');
        beta1_all(:,sub) = beta(1,:);
    end
end
beta1_mean = mean(beta1_all,2,'omitnan');

save('gsr_BOLD_beta_trend_CB.mat', ...
    'beta1_all', 'beta1_mean', 'yeo_all', 'yeo_mean', 'Yeo17_label', '-v7.3');

%% Yeo, beta correlation
clear, close all;
user_addpath(false, false);
addpath('/moonbrainlab/YJPark/matlab/toolbox/ENIGMA-master/matlab/shared/surfaces');

load('gsr_BOLD_beta_trend_CB.mat');

beta1_all_TR = moving_time_window_v2(beta1_all, 21, 105, @(x) x);

yeo_beta_corr_mat = nan([17,22]);
yeo_beta_pval_mat = nan([17,22]);
for sub = 1:22
    try
        [yeo_beta_corr_mat(:,sub), yeo_beta_pval_mat(:,sub)] ...
            = corr(beta1_all_TR(:,sub), yeo_all(:,2:94,sub)');
    end
end
yeo_beta_pval_mat = -log10(yeo_beta_pval_mat);

%%
fig = figure(1);
fig.Position = [100 100 600 500];
for sub = 1:22
    try
        clf;
        hold on; grid on;
        bar(1:17, yeo_beta_corr_mat(:,sub), 'CData',yeo_beta_pval_mat(:,sub), ...
            'FaceColor','flat');
        title(sprintf('sub-%02d',sub));
        cb = colorbar;
        caxis([0 3]);
        cb.Ticks = 0:3;
        cb.TickLabels = 10.^-(0:3);
        cb.Label.String = 'p-value';
        cb.Label.FontSize = 14;
        colormap(flip(slanCM('bone')));
        xticks(1:17);
        xticklabels(Yeo17_label);
        ylim([-0.4 0.4]);

        exportgraphics(fig, sprintf('../figures/checkerboard/gsr_yeo_beta_corr_map_sub-%02d.png',sub));
    end
end

clf;
hold on; grid on;
y = mean(yeo_beta_corr_mat,2,'omitnan');
y_se = std(yeo_beta_corr_mat,0,2,'omitnan')./sqrt(sum(~isnan(yeo_beta_corr_mat),2));
c = mean(yeo_beta_pval_mat,2,'omitnan');
bar(1:17, y, 'CData',c, ...
    'FaceColor','flat');
errorbar(y,y_se, 'LineStyle','none');
title('sub-all');
cb = colorbar;
caxis([0 3]);
cb.Ticks = 0:3;
cb.TickLabels = 10.^-(0:3);
cb.Label.String = 'p-value';
cb.Label.FontSize = 14;
colormap(flip(slanCM('bone')));
xticks(1:17);
xticklabels(Yeo17_label);
ylim([-0.151 0.151]);

exportgraphics(fig, '../figures/checkerboard/gsr_yeo_beta_corr_map_sub-all_v1.png');

clf;
hold on; grid on;
yline(0);
boxplot(yeo_beta_corr_mat', 'PlotStyle','compact', 'Colors','g');
plot(1:17, mean(yeo_beta_corr_mat,2,'omitnan'), 'd', 'MarkerSize',3, ...
    'MarkerFaceColor','b', 'MarkerEdgeColor','k');
box off;
title('sub-all');
xticks(1:17);
xticklabels(Yeo17_label);
ylim([-0.55 0.55]);
legend({'median', 'mean'}, 'NumColumns',2, 'FontSize',14);

exportgraphics(fig, '../figures/checkerboard/gsr_yeo_beta_corr_map_sub-all_v2.png');

%% functional connectivity (gsr, nogsr) 2024.09.10
clear, close all;
user_addpath(false, false);

Yeo = niftiread('../fmri_references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
n_Yeo = 17;
Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};

load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';
file_nogsr = 'func_pp_filter_sm0.mni152.3mm.nii.gz';
file_gsr = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';

sub = 1;
V_nogsr = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
    load_path, sub, sub, file_nogsr));
V_nogsr = zscore(V_nogsr,0,4); V_nogsr = permute(V_nogsr,[4,1,2,3]);

V_gsr = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
    load_path, sub, sub, file_gsr));
V_gsr = zscore(V_gsr,0,4); V_gsr = permute(V_gsr,[4,1,2,3]);

fc_mat_nogsr = nan(n_Yeo);
fc_mat_gsr = nan(n_Yeo);
for y1 = 1:n_Yeo
    idx_tmp1 = find(Yeo == y1);
    for y2 = 1:n_Yeo
        idx_tmp2 = find(Yeo == y2);
        fc_mat_nogsr(y1,y2) = corr(mean(V_nogsr(:,idx_tmp1),2,'omitnan'), ...
            mean(V_nogsr(:,idx_tmp2),2,'omitnan'));
        fc_mat_gsr(y1,y2) = corr(mean(V_gsr(:,idx_tmp1),2,'omitnan'), ...
            mean(V_gsr(:,idx_tmp2),2,'omitnan'));
    end
end

%%
fig = figure(1);
fig.Position = [100 100 1200 400];
clf;
sgtitle(sprintf('sub-%02d', sub), ...
    'FontWeight','bold', 'FontSize',15);
subplot(1,2,1);
heatmap(fc_mat_nogsr, 'YData',Yeo17_label);
caxis([-1 1]);
colormap(flip(slanCM('RdBu')));
title('without GSR');

subplot(1,2,2);
heatmap(fc_mat_gsr, 'YData',Yeo17_label);
caxis([-1 1]);
colormap('jet');
colormap(flip(slanCM('RdBu')));
title('with GSR');
    
% exportgraphics(fig, ...
%     sprintf('../figures/data_check_fig/fmri_issues/comparison_w_gsr_nogsr_sub-%02d.png',sub));

%% functional connectivity (gsr, nogsr) 2024.09.10 YEO17
clear, close all;
user_addpath(false, false);

Yeo = niftiread('../fmri_references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
n_Yeo = 17;
Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};

load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';

fig = figure(1);
fig.Position = [100 100 1300 650];
clf;
% file_name = 'func_pp_filter_sm0.mni152.3mm.nii.gz';
% text(0.9,0.05,'without GSR', 'Units','normalized', ...
%     'FontWeight','bold', 'FontSize',15);
file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';
text(0.9,0.05,'with GSR', 'Units','normalized', ...
    'FontWeight','bold', 'FontSize',15);
axis off;

for y = 1:n_Yeo
    text(0.72,0.18-0.016*y, Yeo17_label{y}, ...
        'HorizontalAlignment','left', 'FontSize',7);
end
% sub = 1;
for sub = 1:22
    try
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
            load_path, sub, sub, file_name));
        V = zscore(V,0,4); V = permute(V,[4,1,2,3]);
        
        fc_mat = nan(n_Yeo);
        for y1 = 1:n_Yeo
            idx_tmp1 = find(Yeo == y1);
            for y2 = 1:n_Yeo
                idx_tmp2 = find(Yeo == y2);
                fc_mat(y1,y2) = corr(mean(V(:,idx_tmp1),2,'omitnan'), ...
                    mean(V(:,idx_tmp2),2,'omitnan'));
            end
        end
        
        subplot_tight(4,6,sub);
        heatmap(fc_mat);
        caxis([-1 1]);
        colormap(flip(slanCM('RdBu')));
        colorbar off;
        title(sprintf('sub-%02d',sub));
    end
end
colorbar;
% exportgraphics(fig, ...
%     '../figures/data_check_fig/fmri_issues/comparison_w_nogsr_YEO17.png');
exportgraphics(fig, ...
    '../figures/data_check_fig/fmri_issues/comparison_w_gsr_YEO17.png');

%% functional connectivity (gsr, nogsr) 2024.09.10 YEO7
clear, close all;
user_addpath(false, false);

Yeo = niftiread('../fmri_references/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
n_Yeo = 7;
Yeo17_label = {'N1: Visual(VIN)', 'N2: Somatomotor(SMN)', 'N3: Dorsal Attention(DAN)', ...
    'N4: Salience(SAN)', 'N5: Limbic(LIN)', 'N6: Frontoparietal(FPN)', ...
    'N7: Default mode network(DMN)'};

load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';

fig = figure(1);
fig.Position = [100 100 1300 650];
clf;
file_name = 'func_pp_filter_sm0.mni152.3mm.nii.gz';
text(0.9,0.05,'without GSR', 'Units','normalized', ...
    'FontWeight','bold', 'FontSize',15);
% file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';
% text(0.9,0.05,'with GSR', 'Units','normalized', ...
%     'FontWeight','bold', 'FontSize',15);
axis off;

for y = 1:n_Yeo
    text(0.72,0.12-0.016*y, Yeo17_label{y}, ...
        'HorizontalAlignment','left', 'FontSize',7);
end
% sub = 1;
for sub = 1:22
    try
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
            load_path, sub, sub, file_name));
        V = zscore(V,0,4); V = permute(V,[4,1,2,3]);
        
        fc_mat = nan(n_Yeo);
        for y1 = 1:n_Yeo
            idx_tmp1 = find(Yeo == y1);
            for y2 = 1:n_Yeo
                idx_tmp2 = find(Yeo == y2);
                fc_mat(y1,y2) = corr(mean(V(:,idx_tmp1),2,'omitnan'), ...
                    mean(V(:,idx_tmp2),2,'omitnan'));
            end
        end
        
        subplot_tight(4,6,sub);
        heatmap(fc_mat);
        caxis([-1 1]);
        colormap(flip(slanCM('RdBu')));
        colorbar off;
        title(sprintf('sub-%02d',sub));
    end
end
colorbar;
exportgraphics(fig, ...
    '../figures/data_check_fig/fmri_issues/comparison_w_nogsr_YEO7.png');
% exportgraphics(fig, ...
%     '../figures/data_check_fig/fmri_issues/comparison_w_gsr_YEO7.png');

%% inter-subject correlation for each network
clear, close all;
user_addpath(false, false);

Yeo = niftiread('../fmri_references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
n_Yeo = 17;
Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};

load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';
n_sub = 22;

file_name = 'func_pp_filter_sm0.mni152.3mm.nii.gz';
title_name = 'without GSR';
% file_name = 'func_pp_filter_gsr_sm0.mni152.3mm.nii.gz';
% title_name = 'with GSR';

inter_subj_mat = nan([n_sub,n_sub,n_Yeo]);
% sub = 1;
for sub1 = 1:n_sub
    try
        V1 = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
            load_path, sub1, sub1, file_name));
        V1 = zscore(V1,0,4); V1 = permute(V1,[4,1,2,3]);

        for sub2 = 1:n_sub
            try
                V2 = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-checker_bold/func_preproc/%s', ...
                    load_path, sub2, sub2, file_name));
                V2 = zscore(V2,0,4); V2 = permute(V2,[4,1,2,3]);


                for y = 1:n_Yeo
                    idx_tmp = find(Yeo == y);
                    inter_subj_mat(sub1,sub2,y) = corr( ...
                        mean(V1(:,idx_tmp),2,'omitnan'), ...
                        mean(V2(:,idx_tmp),2,'omitnan'));
                end    
            end
        end
    end
end


fig = figure(1);
fig.Position = [100 100 1300 700];
clf;
text(0.95,0.05,title_name, 'Units','normalized', ...
    'FontWeight','bold', 'FontSize',15);
axis off;
for y = 1:n_Yeo
    subplot_tight(3,6,y);
    heatmap(inter_subj_mat(:,:,y));
    colorbar;
    colormap(flip(slanCM('RdBu')));
    caxis([-0.5 0.5]);
    title(Yeo17_label{y});
end

exportgraphics(fig, ...
    '../figures/data_check_fig/fmri_issues/inter-subj_corr_w_nogsr_YEO17.png');
% exportgraphics(fig, ...
%     '../figures/data_check_fig/fmri_issues/inter-subj_corr_w_gsr_YEO17.png');

%% inter-subject correlation of beta
clear, close all;

load('yeo_all_and_beta_CB.mat');

corr_mat = nan(n_sub);
for i = 1:n_sub
    for j = 1:n_sub
        corr_mat(i,j) = corr(beta_all_TR(1,:,i)', beta_all_TR(1,:,j)');
    end
end

fig = figure(1);
fig.Position = [100 100 800 600];
clf;
heatmap(corr_mat);
colormap(flip(slanCM('RdBu')));
caxis([-0.5 0.5]);
title('inter-subject corr. of beta');
exportgraphics(fig, ...
    '../figures/data_check_fig/fmri_issues/inter-subj_corr_of_beta.png');