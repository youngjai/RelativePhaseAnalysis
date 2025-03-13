%% Let's try to calculate the relative phase on whole time data
clear, close all;
user_addpath(false,false);

load('/moonbrainlab/raw_data/McDonnell_data/UofM_data/matdata500_notch60_detrand0p5to100/UM_8_detrand0p5to100.mat');

%%
tic;
rel_p = cal_rel_phase_v3(double(data));
toc;


%%
% E11:Fz, E62:Pz, E75:Oz
ch_name = 'E62';
ch_idx = find(matches({chanlocs.labels},ch_name));

fig = figure(1);
fig.Position = [80 100 1800 800];
clf;
hold on; grid on;
plot((1:size(data,1))/fs/60, data(:,ch_idx));
xline(etime/fs/60);
xlim([0,size(data,1)/fs/60]);
ylim([-400 400]);
text(etime/fs/60,300*ones([53,1]),ename, 'Rotation',45);

%%
fig = figure(2);
fig.Position = [200 200 900 700];
clf;
hold on; grid on;
pspectrum(data(:,ch_idx),fs,'power','FrequencyLimits',[0.5 70]);

%%
fig = figure(3);
fig.Position = [80 100 1800 800];
clf;
hold on; grid on;
plot((1:size(rel_p,1))/fs/60, rel_p(:,ch_idx));
xline(etime/fs/60);
xlim([0,size(rel_p,1)/fs/60]);
ylim([-2 2]);
text(etime/fs/60,1.5*ones([53,1]),ename, 'Rotation',45);

