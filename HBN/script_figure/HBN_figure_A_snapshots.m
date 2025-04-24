%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load relaive phase files
% 2. Make temporal snapshots for each group

% 2024. 10. 06. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% relative phase change pattern check 
% ADHD EC
load('D:\Dropbox\MoonBrainLab Raw Data (1)\Healthy Brain Network data\regress_result\Ain_ec_cal_sub.mat')
IDX_ADHD = IDX_sub_Ain_ec;


% CONTROL EC
load('D:\Dropbox\MoonBrainLab Raw Data (1)\Healthy Brain Network data\regress_result\C_ec_new_cal_sub.mat')
IDX_control = IDX_sub;

clear beta* *prop* IDX_sub*



%% load file 
% ADHD EC
% subject1: NDARAY461TZZ_bc_re_EC1_tm50_tw50_sm20_rp 1 session: 8~18 mode4
load('D:\HealthyBrainNetwork_above11\EEGDownload4\Ain_w_124_bp_bc_alpha_ec_original_new\NDARAY461TZZ_bc_re_EC1_tm50_tw50_sm20_rp.mat')
ADHD_topo = topo;

% subject2: NDARAL828WXM_bc_re_EC1_tm50_tw50_sm20_rp 1 session: 70~80 
load('D:\HealthyBrainNetwork_above11\EEGDownload4\C_w_124_bp_bc_alpha_ec_original_new\NDARAM277WZT_bc_re_EC1_tm50_tw50_sm20_rp.mat')
Control_topo = topo;

H = double(H);

% making borderCoords
[dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test(H(1, :), chan_coord_xy(:, 1:2),'smooth',smooth, 'scatter', 1);

%figure condition

formatSpec = '%.2f';

T = 11; % 10 snapshots
interval = 1; % interval : 0.2s
unit = 0.1; % 0.1s

state1 = 'Control';
state2 = 'ADHD';

% % check topo 
% for i = 1:2:length(topo)
%     fig = figure(1)
%     topoplot_figure(Control_topo{i}, borderCoords, xx, yy, Coord, 'scatter', 0, 'line', 0);
%     colormap(slanCM('jet'));
%     clim([-0.75 0.75]);
%     name_con = ['topoplot_control_' num2str(i) '.png'];
%     path_con = ['Y:\etc\contro3_subject3\' name_con]
%     exportgraphics(fig, path_con);
%     clf;
%     
% %     fig = figure(2)
% %     topoplot_figure(ADHD_topo{i}, borderCoords, xx, yy, Coord, 'scatter', 0, 'line', 0)
% %     colormap(slanCM('jet'));
% % %     clim([-0.75 0.75]);
% %     name_adhd = ['topoplot_adhd_' num2str(i) '.png'];
% %     path_adhd = ['Y:\etc\adhd2\' name_con]
% %     exportgraphics(fig, path_adhd);
% %     clf;
% 
% end

start1 = 289; % 46
start2 = 13;

j = 1;
for i  = start2 : 2: (start2+21) 
    topo_set_control{j} = Control_topo{i, 1};
    topo_set_adhd{j} = ADHD_topo{i, 1};
    j = j+1;
end


fig = figure(1);
fig.Position = [50 100 1500 400];

clf;
hold on;

for t = 1:T
    subplot(2,T,t)
    topoplot_figure(topo_set_control{t}, borderCoords, xx, yy, Coord, 'scatter', 0, 'line', 0);
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
    topoplot_figure(topo_set_adhd{t}, borderCoords, xx, yy, Coord, 'scatter', 0, 'line', 0)
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
text(-14.3,2.45,sprintf('Control\n(Eyes-closed)'), 'FontSize',14, ...
    'Units','normalized', 'Rotation',90, 'HorizontalAlignment','center');
text(-14.3,0.45,sprintf('ADHD\n(Eyes-closed)'), 'FontSize',14, ...
    'Units','normalized', 'Rotation',90, 'HorizontalAlignment','center');
cb = colorbar;
cb.Location = 'north';
cb.Position = [0.758 0.444 0.1 0.025];
cb.Ticks = [-0.5 0.5];
cb.TickLabels = {'Lag', 'Lead'};
cb.FontSize = 10;
cb.Label.String = 'Relative phase';
cb.Label.FontSize = 12;
axes(fig, 'Position',[0.753 0.43 0.11 0.15], 'Units','normalized');
axis off;
patch([0 0 1 1], [0 1 1 0], 'w', 'FaceColor','none', ...
    'EdgeColor','k');
xlim([0 1]);
ylim([0 1]);
exportgraphics(fig, 'a_snapshot_of_topoplot.png');
exportgraphics(fig, 'a_snapshot_of_topoplot.eps');
savefig(fig, 'a_snapshot_of_topoplot.fig');