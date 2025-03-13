%% issue 1. average BOLD over time for each voxel is zero
clear, close all;
user_addpath(false,false);

task = 'rest';
load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';
n_sub = 22;

nT = 288; nX = 61*73*61;

V_all = nan([nT,nX,n_sub,2]);
for sub = 1:n_sub
    try
        gsr = '';
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter%s_sm0.mni152.3mm.nii.gz', ...
            load_path,sub,sub,task,gsr));
        V = permute(V, [4,1,2,3]);
        V_all(:,:,sub,1) = V(:,:);

        gsr = '_gsr';
        V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter%s_sm0.mni152.3mm.nii.gz', ...
            load_path,sub,sub,task,gsr));
        V = permute(V, [4,1,2,3]);
        V_all(:,:,sub,2) = V(:,:);
    end
end

fig = figure(1);
fig.Position = [100 100 1600 600];
clf;
subplot(2,1,1);
hold on; grid on;
gsr = '';
title(sprintf('func_pp_filter%s_sm0.mni152.3mm.nii.gz',gsr), ...
    'Interpreter','none');
plot(reshape(mean(V_all(:,:,:,1),1,'omitnan'),1,[]), '.');
xline(1:nX:nX*(n_sub+1));
set(gca, 'FontSize',15);
ylim([-4, 4]*1e-5);
ylabel('$\langle \textsf{BOLD}_i \rangle_T$', ...
    'Interpreter','latex', 'FontSize',18);
xticks(0.5*nX:nX:nX*n_sub);
xticklabels(compose('sub-%02d',1:n_sub));

subplot(2,1,2);
hold on; grid on;
gsr = '_gsr';
title(sprintf('func_pp_filter%s_sm0.mni152.3mm.nii.gz',gsr), ...
    'Interpreter','none');
plot(reshape(mean(V_all(:,:,:,2),1,'omitnan'),1,[]), '.');
xline(1:nX:nX*(n_sub+1));
set(gca, 'FontSize',15);
ylim([-4, 4]*1e-5);
ylabel('$\langle \textsf{BOLD}_i \rangle_T$', ...
    'Interpreter','latex', 'FontSize',18);
xticks(0.5*nX:nX:nX*n_sub);
xticklabels(compose('sub-%02d',1:n_sub));

exportgraphics(fig, '../figures/data_check_fig/fmri_issues/averageBOLD.png');

%% correlation matrix between YEO17 netowrks
clear, close all;
user_addpath(false,false);

n_sub = 22;
task = 'rest';
load_path = '/moonbrainlab/raw_data/NKI_eeg_fmri/preproc';

Yeo = niftiread('../data/references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
Yeo_all = nan([17,288,n_sub]);
Yeo_all_gsr = nan([17,288,n_sub]);
for sub = 1:n_sub
    V_gsr = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_gsr_sm0.mni152.3mm.nii.gz', ...
        load_path,sub,sub,task));
    V_gsr = permute(V_gsr, [4,1,2,3]);
    V = niftiread(sprintf('%s/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_sm0.mni152.3mm.nii.gz', ...
        load_path,sub,sub,task));
    V = permute(V, [4,1,2,3]);
    nT = size(V,1);

    V_arr = nan([17,nT]);
    V_gsr_arr = nan([17,nT]);
    for yi = 1:17
        idx = find(Yeo == yi);
        V_arr(yi,:) = mean(V(:,idx),2);
        V_gsr_arr(yi,:) = mean(V_gsr(:,idx),2);
    end
%     V_arr_z = zscore(V_arr,0,2);
%     V_gsr_arr_z = zscore(V_gsr_arr,0,2);
%     Yeo_all(:,:,sub) = V_arr_z;
%     Yeo_all_gsr(:,:,sub) = V_gsr_arr_z;

    Yeo_all(:,:,sub) = V_arr;
    Yeo_all_gsr(:,:,sub) = V_gsr_arr;
end

Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};

corr_mat = nan([17,17,22,2]);
for sub = 1:n_sub
    corr_mat(:,:,sub,1) = corr(Yeo_all(:,:,sub)',Yeo_all(:,:,sub)');
    corr_mat(:,:,sub,2) = corr(Yeo_all_gsr(:,:,sub)',Yeo_all_gsr(:,:,sub)');
end

%%
fig = figure(1);
fig.Position = [100 100 1200 500];
clf;
subplot(1,2,1);
hold on; grid on;
gsr = '';
title(sprintf('func_pp_filter%s_sm0.mni152.3mm.nii.gz',gsr), ...
    'Interpreter','none', 'FontSize',14);
imagesc(mean(corr_mat(:,:,:,1),3,'omitnan'));
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
xticks(1:17);
xticklabels(compose('N%d',1:17));
yticks(1:17);
yticklabels(Yeo17_label);
caxis([-1 1]);

subplot(1,2,2);
hold on; grid on;
gsr = '_gsr';
title(sprintf('func_pp_filter%s_sm0.mni152.3mm.nii.gz',gsr), ...
    'Interpreter','none', 'FontSize',14);
imagesc(mean(corr_mat(:,:,:,2),3,'omitnan'));
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
xticks(1:17);
xticklabels(compose('N%d',1:17));
yticks(1:17);
yticklabels(Yeo17_label);
caxis([-1 1]);

colormap(flipud(slanCM('Spectral')));
ax = axes(fig, 'Position',[0.77 0.12 0.2 0.8], 'Unit','normalized');
axis off;
cb = colorbar(ax);
cb.Label.String = 'Correlation Coefficient';
cb.Label.FontSize = 15;
cb.Ticks = 0:0.25:1;
cb.TickLabels = -1:0.5:1;

exportgraphics(fig, '../figures/data_check_fig/fmri_issues/yeo17_correlation_matrix.png');

%%
clear, close all;

sub = 5;
task = 'rest';
V_gsr = niftiread(sprintf('/moonbrainlab/raw_data/NKI_eeg_fmri/preproc/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_gsr_sm0.mni152.3mm.nii.gz',sub,sub,task));
V_gsr = permute(V_gsr, [4,1,2,3]);
V = niftiread(sprintf('/moonbrainlab/raw_data/NKI_eeg_fmri/preproc/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_sm0.mni152.3mm.nii.gz',sub,sub,task));
V = permute(V, [4,1,2,3]);
nT = size(V,1);

fig = figure(1);
fig.Position = [100 100 1800 400];
clf;
hold on; grid on;
plot(V(:,30,30,30));
plot(V_gsr(:,30,30,30));
legend({'sm0', 'gsr_sm0'}, 'Interpreter','none');

mean(abs(V-V_gsr), 'all', 'omitnan')


%%
Yeo = niftiread('../data/references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
V_arr = nan([17,nT]);
V_gsr_arr = nan([17,nT]);
for yi = 1:17
    idx = find(Yeo == yi);
    V_arr(yi,:) = mean(V(:,idx),2);
    V_gsr_arr(yi,:) = mean(V_gsr(:,idx),2);
end
V_arr_z = zscore(V_arr,0,2);
V_gsr_arr_z = zscore(V_gsr_arr,0,2);

corr_mat = corr(V_arr_z',V_arr_z');        
corr_mat_gsr = corr(V_gsr_arr_z',V_gsr_arr_z');        

%%
fig2 = figure(21);
fig2.Position = [100 100 1500 600];
clf;
subplot(1,2,1);
hold on; grid on;
title('sm0', 'Interpreter','none');
imagesc(corr_mat);
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
colorbar;
caxis([-1, 1]);

subplot(1,2,2);
hold on; grid on;
title('gsr_sm0', 'Interpreter','none');
imagesc(corr_mat_gsr);
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
colorbar;
caxis([-1, 1]);

%%
clear, close all;

n_sub = 22;
% task = 'rest';
task = 'checker';

Yeo = niftiread('../data/references/Yeo2011_17Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');
% Yeo_all = nan([17,288,n_sub]);
% Yeo_all_gsr = nan([17,288,n_sub]);
Yeo_all = nan([17,98,n_sub]);
Yeo_all_gsr = nan([17,98,n_sub]);
for sub = 1:n_sub
    V_gsr = niftiread(sprintf('/moonbrainlab/raw_data/NKI_eeg_fmri/preproc/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_gsr_sm0.mni152.3mm.nii.gz',sub,sub,task));
    V_gsr = permute(V_gsr, [4,1,2,3]);
    V = niftiread(sprintf('/moonbrainlab/raw_data/NKI_eeg_fmri/preproc/sub-%02d/ses-01/func/sub-%02d_ses-01_task-%s_bold/func_preproc/func_pp_filter_sm0.mni152.3mm.nii.gz',sub,sub,task));
    V = permute(V, [4,1,2,3]);
    nT = size(V,1);

    V_arr = nan([17,nT]);
    V_gsr_arr = nan([17,nT]);
    for yi = 1:17
        idx = find(Yeo == yi);
        V_arr(yi,:) = mean(V(:,idx),2);
        V_gsr_arr(yi,:) = mean(V_gsr(:,idx),2);
    end
    V_arr_z = zscore(V_arr,0,2);
    V_gsr_arr_z = zscore(V_gsr_arr,0,2);
%     V_arr_z = V_arr./std(V_arr,0,2);
%     V_gsr_arr_z = V_gsr_arr./std(V_arr,0,2);

    Yeo_all(:,:,sub) = V_arr_z;
    Yeo_all_gsr(:,:,sub) = V_gsr_arr_z;
end
save(sprintf('Yeo_all_%s.mat',task), ...
    'n_sub', 'Yeo', 'Yeo_all', 'Yeo_all_gsr', '-v7.3');

%%
clear, close all;

task = 'checker';
load(sprintf('Yeo_all_%s.mat',task));
corr_mat = nan([17,17,n_sub]);
corr_mat_gsr = nan([17,17,n_sub]);
for sub = 1:n_sub
    corr_mat(:,:,sub) = corr(Yeo_all(:,:,sub)',Yeo_all(:,:,sub)');        
    corr_mat_gsr(:,:,sub) = corr(Yeo_all_gsr(:,:,sub)',Yeo_all_gsr(:,:,sub)');    
end

%%
sub = 1;
yeo_idx1 = 1;
yeo_idx2 = 15;
figure(23);
clf;
subplot(2,1,1);
hold on; grid on;
plot((1:288)*2.1, Yeo_all(yeo_idx1,:,sub));
plot((1:288)*2.1, Yeo_all(yeo_idx2,:,sub));
legend({'VisCent', 'DefaultB'}, ...
    'Interpreter','none');
mean(Yeo_all(yeo_idx1,:,sub))
corr(Yeo_all(yeo_idx1,:,sub)',Yeo_all(yeo_idx2,:,sub)')

subplot(2,1,2);
hold on; grid on;
plot((1:288)*2.1, Yeo_all_gsr(yeo_idx1,:,sub));
plot((1:288)*2.1, Yeo_all_gsr(yeo_idx2,:,sub));
legend({'VisCent', 'DefaultB'}, ...
    'Interpreter','none');
mean(Yeo_all_gsr(yeo_idx1,:,sub))
corr(Yeo_all_gsr(yeo_idx1,:,sub)',Yeo_all_gsr(yeo_idx2,:,sub)')

%%
figure(243);
clf;
subplot(1,2,1);
hold on; grid on;
histogram(Yeo_all(yeo_idx1,:,sub),-3:0.1:3, 'Normalization','pdf');

subplot(1,2,2);
hold on; grid on;
histogram(Yeo_all_gsr(yeo_idx1,:,sub),-3:0.1:3, 'Normalization','pdf');

%%
sub = 12;

fig2 = figure(21);
fig2.Position = [100 100 1500 600];
clf;
subplot(1,2,1);
hold on; grid on;
title('sm0', 'Interpreter','none');
imagesc(mean(corr_mat,3));
% corr_mat_z = zscore(corr_mat,0,[1,2]);
% imagesc(mean(corr_mat_z,3));
% imagesc(corr_mat(:,:,sub));
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
colorbar;
% caxis([-1, 1]);
caxis([0, 1]);

subplot(1,2,2);
hold on; grid on;
title('gsr_sm0', 'Interpreter','none');
imagesc(mean(corr_mat_gsr,3));
% corr_mat_gsr_z = zscore(corr_mat_gsr,0,[1,2]);
% imagesc(mean(corr_mat_gsr_z,3));
% imagesc(corr_mat_gsr(:,:,sub));
set(gca, 'YDir','reverse');
xlim([0.5 17.5]);
ylim([0.5 17.5]);
colorbar;
caxis([-1, 1]);

Yeo17_label = {'N1: VisCent', 'N2: VisPeri', 'N3: SomMotA', 'N4: SomMotB', ...
    'N5: DorsAttnA', 'N6: DorsAttnB', 'N7: SalVentAttnA', 'N8: SalVentAttnB', ...
    'N9: LimbicB', 'N10: LimbicA', 'N11: ContA', 'N12: ContB', ...
    'N13: ContC', 'N14: DefaultA', 'N15: DefaultB', 'N16: DefaultC', ...
    'N17: TempPar'};

%%
load('../data/references/all_yeo_and_beta_EO.mat','beta_all_TR');

figure(1412);
clf;
subplot(2,1,1);
hold on; grid on;
y1 = reshape(Yeo_all(1,2:284,:), 1,[]);
y2 = reshape(Yeo_all(15,2:284,:), 1,[]);
b1 = reshape(zscore(beta_all_TR(1,:,:),0,2), 1,[]);
[~,idx_sort] = sort(b1);
plot(smooth(y1(idx_sort),100));
plot(smooth(y2(idx_sort),100));
plot(smooth(b1(idx_sort),100));

subplot(2,1,2);
hold on; grid on;
y1 = reshape(Yeo_all_gsr(1,2:284,:), 1,[]);
y2 = reshape(Yeo_all_gsr(15,2:284,:), 1,[]);
[~,idx_sort] = sort(b1);
plot(smooth(y1(idx_sort),100));
plot(smooth(y2(idx_sort),100));
plot(smooth(b1(idx_sort),100));


%%
load('../data/references/all_yeo_and_beta_EO.mat','beta_all_TR');

corr_mat = nan([17,22]);
for sub = 1:22
    corr_mat(:,sub) = corr(squeeze(beta_all_TR(1,:,sub))',squeeze(Yeo_all(:,1:283,sub))');
end

corr_mat_z = zscore(corr_mat,0,1);
figure(4235);
clf;
subplot(1,3,1);
hold on; grid on;
bar(mean(corr_mat_z,2));

subplot(1,3,2);
hold on; grid on;
bar(mean(corr_mat,2));
% ylim([-0.2 0.2]);

figure(124143);
clf;
for sub = 1:22
    subplot(5,5,sub);
    hold on; grid on;
    bar(corr_mat(:,sub));
end

%%
yeo_trend = zeros([n_sub,283]);
beta_trend = zeros([n_sub,283]);

for sub = 1:n_sub
    [~,idx_s] = sort(Yeo_all(15,1:283,sub));
    yeo_trend(sub,:) = Yeo_all(15,idx_s,sub);
    beta_trend(sub,:) = beta_all_TR(2,idx_s,sub);
end

%% trend
beta_trend = zscore(beta_trend,0,2);
sub = 20;
figure(15343);
clf;
hold on; grid on;
    
plot((1:283)*2.1,yeo_trend(sub,:));
plot((1:283)*2.1,beta_trend(sub,:));

% plot((1:283)*2.1,mean(yeo_trend,1));
% plot((1:283)*2.1,mean(beta_trend,1));
