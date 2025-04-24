%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load relative phase files (weighted)
% 2. Perform statistical analysis (t-test, FDR)
% 3. Create figures for mode and dwell time distributions
% 4. Plot results with violin plots and statistical significance (brackets and stars)
% 
% 2025. 03. 12. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% 100 ms weighted

% load the files
file_path = 'D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result\';
load([file_path 'regress_result_new_cal_sub']);

condition = 'alpha';  % 'theta' or 'delta'

files = {
    ['C_ec'],
    ['C_eo'],
    ['Ain_ec'],
    ['Ain_eo']
};

% load the files
for i = 1:length(files)
    file_name = files{i};
    file_task = file_name;

    var_name = ['w_prop_sub_' file_task ];       


    eval([' w_prop_sub =' var_name ';']);  
    
    for subj_no = 1:length(w_prop_sub)
        
        data = w_prop_sub(subj_no, 1:5);
    
        for j = 1:5
        
            mode(j, 1:4) = [data(j).avg_weight];
            dwell(j, 1:4) = [data(j).dwell_time];
    
        end
    
        sub_mode(subj_no, 1:4) = mean(mode);
        sub_dwell(subj_no, 1:4) = mean(dwell)*0.01;
    
    end

    var_name_mode = ['sub_mode_' file_task]; 
    var_name_dwell = ['sub_dwell_' file_task]; 

    eval([var_name_mode ' = sub_mode;']);  
    eval([var_name_dwell ' = sub_dwell;']);  

    mean_mode = mean(sub_mode); % group
    mean_dwell = mean(sub_dwell); % group

    se_mode = std(sub_mode)/sqrt(length(sub_mode));
    se_dwell = std(sub_dwell)/sqrt(length(sub_dwell));

    var_name_mode_m = ['mean_mode_' file_task]; 
    var_name_dwell_m = ['mean_dwell_' file_task]; 

    eval([var_name_mode_m ' = mean_mode;']);  
    eval([var_name_dwell_m ' = mean_dwell;']);  

    var_name_mode_se = ['se_mode_' file_task]; 
    var_name_dwell_se = ['se_dwell_' file_task]; 

    eval([var_name_mode_se ' = se_mode;']);  
    eval([var_name_dwell_se ' = se_dwell;']);  

    clear sub_mode sub_dwell mode dwell data w_prop_sub mean_mode mean_dwell se_mode se_dwell

end


%% t-test
[h_ec_mode, p_ec_mode, ci_ec_mode, stats_ec_mode] = ttest2(sub_mode_C_ec, sub_mode_Ain_ec);
[h_eo_mode, p_eo_mode, ci_eo_mode, stats_eo_mode] = ttest2(sub_mode_C_eo, sub_mode_Ain_eo);
[h_ec_dwell, p_ec_dwell, ci_ec_dwell, stats_ec_dwell] = ttest2(sub_dwell_C_ec, sub_dwell_Ain_ec);
[h_eo_dwell, p_eo_dwell, ci_eo_dwell, stats_eo_dwell] = ttest2(sub_dwell_C_eo, sub_dwell_Ain_eo);


%% FDR
% p_values = [p_eo, p_ec];
% q_values = mafdr(p_values, 'BHFDR', true);

%% figure 

load('D:\HealthyBrainNetwork_above11\code\newcode\mode_colors.mat')

mode = vertcat(mean_mode_C_eo, mean_mode_Ain_eo, mean_mode_C_ec, mean_mode_Ain_ec);
mode_re = horzcat(mode(1:4, 1)', mode(1:4, 2)', mode(1:4, 3)',  mode(1:4, 4)');

mode_se = vertcat(se_mode_C_eo, se_mode_Ain_eo, se_mode_C_ec, se_mode_Ain_ec);
mode_se_re = horzcat(mode_se(1:4, 1)', mode_se(1:4, 2)', mode_se(1:4, 3)',  mode_se(1:4, 4)');

dwell = vertcat(mean_dwell_C_eo, mean_dwell_Ain_eo, mean_dwell_C_ec, mean_dwell_Ain_ec);
dwell_re = horzcat(dwell(1:4, 1)', dwell(1:4, 2)', dwell(1:4, 3)',  dwell(1:4, 4)');

dwell_se = vertcat(se_dwell_C_eo, se_dwell_Ain_eo, se_dwell_C_ec, se_dwell_Ain_ec);
dwell_se_re = horzcat(dwell_se(1:4, 1)', dwell_se(1:4, 2)', dwell_se(1:4, 3)',  dwell_se(1:4, 4)');


%% figure

for i = 1:length(c_list)

    num = (i-1)*4+1;

    data_dwell{num} = sub_dwell_C_eo(:, i);
    data_dwell{num+1} = sub_dwell_Ain_eo(:, i);
    data_dwell{num+2} = sub_dwell_C_ec(:, i);
    data_dwell{num+3} = sub_dwell_Ain_ec(:, i);

    data_mode{num} = sub_mode_C_eo(:, i);
    data_mode{num+1} = sub_mode_Ain_eo(:, i);
    data_mode{num+2} = sub_mode_C_ec(:, i);
    data_mode{num+3} = sub_mode_Ain_ec(:, i);

end

mean_mode = cellfun(@(x) mean(x, 'all', 'omitnan'), data_mode);
mean_dwell = cellfun(@(x) mean(x, 'all', 'omitnan'), data_dwell);

%%

fig = figure(1);
fig.Position = [100, 100, 1600, 500];
clf;
ax_main = axes(fig, 'Position', [0.13 0.25 0.7750 0.5705]);
hold on; grid on;

hold on; grid on;

% Shading for intervals
shading_intervals = {[3, 6], [10, 13], [17, 20], [24, 27]};
% y_limits = get(gca, 'YLim');  % Get the y-axis limits
y_limits = [0 80];
for i = 1:length(shading_intervals)
    x = shading_intervals{i};
    x_patch = [x(1), x(2), x(2), x(1)];  % x-coordinates for the shading rectangle
    y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];  % y-coordinates for the shading rectangle
    patch(x_patch, y_patch, [1, 0.95, 0.46], 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Add shading with light gray color
end

% Define bar positions and groups
x_positions = [1, 2, 4, 5, 8, 9, 11, 12, 15, 16, 18, 19, 22, 23, 25, 26];
colors = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
alphas = [0.5, 0.6, 0.5, 0.6];

for i = 1:16  
    % Find the alpha value for the respective color
    alpha_idx = mod(i-1, 4) + 1;
    alpha_value = alphas(alpha_idx); 
    
    % Run daviolinplot
    h = daviolinplot(data_mode{i}, 'violin', 'full', ...
        'colors', c_list(colors(i), :), ...
        'violinalpha', alpha_value, ... 
        'boxwidth', 1.2, 'outsymbol', 'rx');

    
    xlim([0 27])
    xticks(1:27);
    xticklabels(string(1:27));
    ylim([0 80])
    
    shift = x_positions(i) - 1; % Move from default position (i) to new position
    h.bx.XData = h.bx.XData + shift;
    h.md.XData = h.md.XData + shift;
    h.ds.XData = h.ds.XData + shift;
    h.ot.XData = h.ot.XData + shift;
    for j = 1:size(h.wh, 3)
        h.wh(1,1,j).XData = h.wh(1,1,j).XData + shift;
    end
end


% Set the X-axis
xticks(x_positions);
set(gca, 'XTick', [1:2, 4:5, 8:9, 11:12, 15:16, 18:19, 22:23, 25:26], 'XTickLabel', {'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD'}, 'FontSize', 16);

ylabel('Normalized Occurrence (%)');



% Dotted lines 

for j = 1:8

    s = (j-1)*2 + 1;
    e = s+1;

    plot(x_positions(s:e), mean_mode(s:e), ':', ...
        'LineWidth', 4, 'Color', c_list(colors(s),:));
   
end

% Add horizontal bracket and stars
x_bracket = [22, 23, 25, 26];  % X-axis positions for the bracket
y_bracket = [50 70];        % Y-axis position for the bracket


line([x_bracket(1), x_bracket(2)], [y_bracket(1), y_bracket(1)], 'Color', 'k', 'LineWidth', 1.5);  
line([x_bracket(1), x_bracket(1)], [y_bracket(1), y_bracket(1)-2], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(2), x_bracket(2)], [y_bracket(1), y_bracket(1)-2], 'Color', 'k', 'LineWidth', 1.5); 

line([x_bracket(3), x_bracket(4)], [y_bracket(2), y_bracket(2)], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(3), x_bracket(3)], [y_bracket(2), y_bracket(2)-2], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(4), x_bracket(4)], [y_bracket(2), y_bracket(2)-2], 'Color', 'k', 'LineWidth', 1.5); 


% Add stars above the bracket
text(mean(x_bracket(1:2)), y_bracket(1) + 3, '***', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
text(mean(x_bracket(3:4)), y_bracket(2) + 3, '***', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');


%% 100ms dwell

fig = figure(2);
fig.Position = [100, 100, 1600, 500];
clf;
ax_main = axes(fig, 'Position', [0.13 0.25 0.7750 0.5705]);
hold on; grid on;

hold on; grid on;
% Shading for intervals
shading_intervals = {[3, 6], [10, 13], [17, 20], [24, 27]};
% y_limits = get(gca, 'YLim');  % Get the y-axis limits
y_limits = [0 1];
for i = 1:length(shading_intervals)
    x = shading_intervals{i};
    x_patch = [x(1), x(2), x(2), x(1)];  % x-coordinates for the shading rectangle
    y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];  % y-coordinates for the shading rectangle
    patch(x_patch, y_patch, [1, 0.95, 0.46], 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Add shading with light gray color
end

% Define bar positions and groups
x_positions = [1, 2, 4, 5, 8, 9, 11, 12, 15, 16, 18, 19, 22, 23, 25, 26];
colors = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
alphas = [0.5, 0.6, 0.5, 0.6];

for i = 1:16  
    % Find the alpha value for the respective color
    alpha_idx = mod(i-1, 4) + 1; 
    alpha_value = alphas(alpha_idx); 
    
    % Run daviolinplot
    h = daviolinplot(data_dwell{i}, 'violin', 'full', ...
        'colors', c_list(colors(i), :), ...
        'violinalpha', alpha_value, ... 
        'boxwidth', 1.2, 'outsymbol', 'rx');


    
    xlim([0 27])
    xticks(1:27);
    xticklabels(string(1:27));
    ylim([0 0.65])
    
    shift = x_positions(i) - 1;
    h.bx.XData = h.bx.XData + shift;
    h.md.XData = h.md.XData + shift;
    h.ds.XData = h.ds.XData + shift;
    h.ot.XData = h.ot.XData + shift;
    for j = 1:size(h.wh, 3)
        h.wh(1,1,j).XData = h.wh(1,1,j).XData + shift;
    end
end


% Set the X-axis
xticks(x_positions);
set(gca, 'XTick', [1:2, 4:5, 8:9, 11:12, 15:16, 18:19, 22:23, 25:26], 'XTickLabel', {'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD'}, 'FontSize', 16);

ylabel('Normalized dwell time (s)');



% Dotted lines 

for j = 1:8

    s = (j-1)*2 + 1;
    e = s+1;

    plot(x_positions(s:e), mean_dwell(s:e), ':', ...
        'LineWidth', 4, 'Color', c_list(colors(s),:));
   
end

% Add horizontal bracket and stars
x_bracket = [22, 23, 25, 26]; 
y_bracket = [0.35 0.57];     


line([x_bracket(1), x_bracket(2)], [y_bracket(1), y_bracket(1)], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(1), x_bracket(1)], [y_bracket(1), y_bracket(1)-0.02], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(2), x_bracket(2)], [y_bracket(1), y_bracket(1)-0.02], 'Color', 'k', 'LineWidth', 1.5); 

line([x_bracket(3), x_bracket(4)], [y_bracket(2), y_bracket(2)], 'Color', 'k', 'LineWidth', 1.5);  
line([x_bracket(3), x_bracket(3)], [y_bracket(2), y_bracket(2)-0.02], 'Color', 'k', 'LineWidth', 1.5); 
line([x_bracket(4), x_bracket(4)], [y_bracket(2), y_bracket(2)-0.02], 'Color', 'k', 'LineWidth', 1.5); 



% Add stars above the bracket
text(mean(x_bracket(1:2)), y_bracket(1) + 0.03, '***', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
text(mean(x_bracket(3:4)), y_bracket(2) + 0.03, '***', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');

