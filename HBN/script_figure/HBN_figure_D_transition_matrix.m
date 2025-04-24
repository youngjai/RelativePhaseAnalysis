
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load relative phase transition probability data
% 2. Create a figure for the transition matrix
% 
% 2025. 02. 27. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd('D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result')
load("regress_result_new_cal_sub.mat")

K = 4; unit = 0.1;

% load transition matrix

% trans_mat_1 == control / 2 == adhd
n_sub_1 = length(prop_sub_C_ec);
n_sub_2 = length(prop_sub_Ain_ec);

trans_mat_1 = nan([K,K,n_sub_1,1]);
trans_mat_2 = nan([K,K,n_sub_2,1]);

for sub = 1:n_sub_1
    try
        % control group
        period_mat_1 = nan([K,K,5,1]);
        for period = 1:5
            period_mat_1(:, :, period, 1) = prop_sub_C_ec(sub, period).trans_prob_mat;
        end
        
        trans_mat_1(:,:,sub,1) = mean(period_mat_1, 3, 'omitnan');

%       Control group

        period_mat_2 = nan([K,K,5,1]);
        if sub < 41
            for period = 1:5
                period_mat_2(:, :, period, 1) = prop_sub_Ain_ec(sub, period).trans_prob_mat;
            end 

            trans_mat_2(:,:,sub,1) = mean(period_mat_2, 3, 'omitnan');

        end
%       ADHD group
    
    end
end


% ec_t_b : ec_transition_beta
% eo_t_b : eo_transition_beta


% EC-ADHD 
   
% calculate mean difference 
mean_control_ec = mean(trans_mat_1, 3, 'omitnan');  % calculate the mean for each column in the control group (ignoring NaN)
mean_ADHD_ec = mean(trans_mat_2, 3, 'omitnan');        % calculate the mean for each column in the ADHD group (ignoring NaN)

% mean difference (control - ADHD)
trans_mat_mean_diff = mean_ADHD_ec - mean_control_ec;

% calculate relative percent change
relative_change = (trans_mat_mean_diff / mean_control_ec) * 100;

% perform independent sample t-test
% initialize array to store results (4x4 size)
ttest2_mat = nan(4, 4);

% perform t-test for each row (i) and column (j)
for i = 1:4
    for j = 1:4
        % extract values for the i, j position from the first group (trans_mat_1)
        val_1 = squeeze(trans_mat_1(i, j, :));
        
        % extract values for the i, j position from the second group (trans_mat_2)
        val_2 = squeeze(trans_mat_2(i, j, :));
        
        % remove NaN values
        val_1 = val_1(~isnan(val_1));
        val_2 = val_2(~isnan(val_2));
        
        % perform t-test and store p-value
        [~, p_value] = ttest2(val_1, val_2);
        ttest2_mat(i, j) = p_value;  % Store p-value
    end
end

%% matrix figure
margin_height = 0.16;
margin_width = 0.04;
fig = figure(1);
fig.Position = [100 100 1600 500];
clf;

% draw heatmap
subplot_tight(1,3,1,[margin_height,margin_width]);
hold on; grid on;

% plot image: -log10 transformed t-test results, value range [1, 3]
imagesc(-log10(ttest2_mat), [1,3]);
set(gca, 'YDir','reverse');
colormap(gca, flip(slanCM('bone')));

cb = colorbar;
cb.Label.String = 'p-value';
cb.Label.FontSize = 18;
cb.Ticks = 1:3;
cb.TickLabels = 10.^-(1:3);
set(gca, 'FontSize', 15);

xlabel('Mode (to)', 'FontSize', 18);
xticks(1:4); 
xlim([0.5, 4.5]);
ylabel('Mode (from)', 'FontSize', 18);
yticks(1:4); 
ylim([0.5, 4.5]);

[xTxt, yTxt] = ndgrid(1:4, 1:4); % 
text(xTxt(:), yTxt(:), num2cell(reshape(round(relative_change*1e4)./1e4,1,[])'), ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'FontSize', 15, 'Color', 'w');

xline(0.5:4.5); % Set xline for 4x4 matrix
yline(0.5:4.5); % Set yline for 4x4 matrix

title('Relative Change (%)', 'FontSize', 15, 'FontWeight', 'normal');

% save path and file name settings
save_path = 'D:\Dropbox\Moonbrainlab\Manuscripts\Relative Phase Project\figures\fig files\fig 4 ADHD';
file_name = 'd_transition_figure_unweighted';

full_file_path = fullfile(save_path, [file_name, '.svg']);

% save current figure as SVG
fig = gcf;  % Get current figure handle
print(fig, full_file_path, '-dsvg');


% figure 2, 3

% initial setup and graph generation
state = 'Control';
trans_mean = mean_control_ec;
trans_mean_re = trans_mean';
G = digraph(trans_mean_re);

% create subplot and adjust position
ax = subplot_tight(1,3,2,[margin_height,margin_width]);
ax.Position(1) = ax.Position(1)+0.04;
ax.Position(2) = ax.Position(2)-0.1;
ax.Position(4) = ax.Position(4)+0.08;
hold on; grid on;

xlim([-1.6 1.6]);
ylim([-1.7 1.6]);
axis off;

% set node positions (4 nodes)
x = [0 -1 1  0]; 
y = [1  0 0 -1];
xwidth = 0.3; 
ywidth = 0.4;
x_img_l = x - xwidth; 
x_img_r = x + xwidth;
y_img_d = y - ywidth; 
y_img_u = y + ywidth;

% Load and display images
modes = cell([1, 4]);
for k = 1:4
    modes{k} = imread(['C' num2str(k) '.png']);
    imagesc([x_img_l(k) x_img_r(k)], [y_img_u(k) y_img_d(k)], modes{k});
end

% First graph plot (no arrows)
p = plot(G);
p.XData = x; 
p.YData = y;
p.Marker = 'none';
p.NodeLabel = repmat({''}, 1, numnodes(G));
p.LineWidth = 10 * G.Edges.Weight;
p.ArrowSize = 0 * G.Edges.Weight;
colormap(gca, flipud(gray)); 
clim([0, 0.4]);
p.EdgeCData = G.Edges.Weight; 
p.EdgeAlpha = 1;

% Create graph again after removing diagonal elements
G = digraph(trans_mean_re - diag(diag(trans_mean_re)));
p = plot(G);
p.XData = x; 
p.YData = y;
p.Marker = '.'; 
p.MarkerSize = 25 * ones([1, 4]); 
p.NodeColor = [60 60 70] / 255;
p.NodeLabel = repmat({''}, 1, numnodes(G));
p.LineWidth = 10 * G.Edges.Weight;
p.ArrowSize = 60 * G.Edges.Weight; 
p.ArrowPosition = 0.65;
colormap(gca, flipud(gray)); 
clim([0, 0.4]);
p.EdgeCData = G.Edges.Weight; 
p.EdgeAlpha = 1;

% set node labels (4 modes)
text(p.XData, p.YData + 0.47, {'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4'}, ...
    'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);

% add state text
text(0, gca().YLim(2) + 0.2, 'Control', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');
% ---------------------------------------------------------------
% ADHD group
state = 'ADHD';
trans_mean2 = mean_ADHD_ec;
trans_mean2_re = trans_mean2';
G = digraph(trans_mean2_re);

% create subplot and adjust position
ax = subplot_tight(1,3,3,[margin_height,margin_width]);
ax.Position(2) = ax.Position(2)-0.1;
ax.Position(4) = ax.Position(4)+0.08;
hold on; grid on;

% set graph range
xlim([-1.6 1.6]);
ylim([-1.7 1.6]);
axis off;

% set node positions (4 nodes)
x = [0 -1 1  0]; 
y = [1  0 0 -1];
xwidth = 0.3; 
ywidth = 0.4;
x_img_l = x - xwidth; 
x_img_r = x + xwidth;
y_img_d = y - ywidth; 
y_img_u = y + ywidth;

% load and display images
modes = cell([1, 4]);
for k = 1:4
    modes{k} = imread(['C' num2str(k) '.png']);
    imagesc([x_img_l(k) x_img_r(k)], [y_img_u(k) y_img_d(k)], modes{k});
end

% first graph plot (no arrows)
p = plot(G);
p.XData = x; 
p.YData = y;
p.Marker = 'none';
p.NodeLabel = repmat({''}, 1, numnodes(G));
p.LineWidth = 10 * G.Edges.Weight;
p.ArrowSize = 0 * G.Edges.Weight;
colormap(gca, flipud(gray)); 
clim([0, 0.4]);
p.EdgeCData = G.Edges.Weight; 
p.EdgeAlpha = 1;

% create graph again after removing diagonal elements
G = digraph(trans_mean2_re - diag(diag(trans_mean2_re)));
p = plot(G);
p.XData = x; 
p.YData = y;
p.Marker = '.'; 
p.MarkerSize = 25 * ones([1, 4]); 
p.NodeColor = [60 60 70] / 255;
p.NodeLabel = repmat({''}, 1, numnodes(G));
p.LineWidth = 10 * G.Edges.Weight;
p.ArrowSize = 60 * G.Edges.Weight; 
p.ArrowPosition = 0.65;
colormap(gca, flipud(gray)); 
clim([0, 0.4]);
p.EdgeCData = G.Edges.Weight; 
p.EdgeAlpha = 1;

% set node labels (4 modes)
text(p.XData, p.YData + 0.47, {'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4'}, ...
    'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);

% add state text
text(0, gca().YLim(2) + 0.2, 'ADHD', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');