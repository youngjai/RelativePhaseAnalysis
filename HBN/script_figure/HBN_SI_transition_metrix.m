
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load relative phase transition probability data
% 2. Visualize results as heatmaps with annotated transition probabilities
% 
% 2025. 03. 16. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 100ms
% transition matrix

cd('D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result')
load("regress_result_new_cal_sub.mat")

K = 4; unit = 0.1;

% Load transition matrix

% trans_mat_1 == control / 2 == adhd
n_sub_1 = length(prop_sub_C_eo);
n_sub_2 = length(prop_sub_Ain_eo);
n_sub_3 = length(prop_sub_C_ec);
n_sub_4 = length(prop_sub_Ain_ec);

trans_mat_1 = nan([K,K,n_sub_1,1]);
trans_mat_2 = nan([K,K,n_sub_2,1]);
trans_mat_3 = nan([K,K,n_sub_3,1]);
trans_mat_4 = nan([K,K,n_sub_4,1]);

for sub = 1:n_sub_1
    try
        % control_eo
        period_mat_1 = nan([K,K,5,1]);
        for period = 1:5
            period_mat_1(:, :, period, 1) = prop_sub_C_eo(sub, period).trans_prob_mat;
        end
        
        trans_mat_1(:,:,sub,1) = mean(period_mat_1, 3, 'omitnan');

        % adhd_eo
        period_mat_2 = nan([K,K,5,1]);
        if sub < 41
            for period = 1:5
                period_mat_2(:, :, period, 1) = prop_sub_Ain_eo(sub, period).trans_prob_mat;
            end 

            trans_mat_2(:,:,sub,1) = mean(period_mat_2, 3, 'omitnan');
        end

        % control_ec
        period_mat_3 = nan([K,K,5,1]);
        for period = 1:5
            period_mat_3(:, :, period, 1) = prop_sub_C_ec(sub, period).trans_prob_mat;
        end
        
        trans_mat_3(:,:,sub,1) = mean(period_mat_3, 3, 'omitnan');

        % adhd_ec
        period_mat_4 = nan([K,K,5,1]);
        if sub < 41
            for period = 1:5
                period_mat_4(:, :, period, 1) = prop_sub_Ain_ec(sub, period).trans_prob_mat;
            end 

            trans_mat_4(:,:,sub,1) = mean(period_mat_4, 3, 'omitnan');
        end
    end
end

% Average transition matrices (excluding NaNs)
mean_control_eo = mean(trans_mat_1, 3, 'omitnan');  
mean_ADHD_eo = mean(trans_mat_2, 3, 'omitnan');      
mean_control_ec = mean(trans_mat_3, 3, 'omitnan');  
mean_ADHD_ec = mean(trans_mat_4, 3, 'omitnan');   

% Perform t-tests for each matrix element (i,j)
for i = 1:4
    for j = 1:4
        val_1 = squeeze(trans_mat_1(i, j, :));
        val_2 = squeeze(trans_mat_2(i, j, :));
        val_3 = squeeze(trans_mat_3(i, j, :));
        val_4 = squeeze(trans_mat_4(i, j, :));
        
        val_1 = val_1(~isnan(val_1));
        val_2 = val_2(~isnan(val_2));
        val_3 = val_3(~isnan(val_3));
        val_4 = val_4(~isnan(val_4));
        
        [h_eo, p_value_eo, ci_eo, stats_eo] = ttest2(val_1, val_2);
        ttest2_p_eo(i, j) = p_value_eo;
        ttest2_h_eo(i, j) = h_eo;
        ttest2_stats_eo(i, j) = stats_eo;

        [h_ec, p_value_ec, ci_ec, stats_ec] = ttest2(val_3, val_4);
        ttest2_p_ec(i, j) = p_value_ec;
        ttest2_h_ec(i, j) = h_ec;  
        ttest2_stats_ec(i, j) = stats_ec;
    end
end

%% Matrix figure setup
margin_height = 0.16;
margin_width = 0.04;
fig = figure(1);
fig.Position = [100 100 2000 500];
clf;

transitions = {'mean_control_eo', 'mean_ADHD_eo', ...
    'mean_control_ec', 'mean_ADHD_ec'};

for i = 1:4
    subplot_tight(1,4,i,[margin_height,margin_width]);
    hold on; grid on;

    data = eval(transitions{i});
    
    % Draw heatmap: display transition probability matrix
    imagesc(data, [0, 1]);
    set(gca, 'YDir','reverse');
    colormap('hot');
    
    caxis([0 1])
    cb = colorbar;
    cb.Ticks = [0 0.5 1];
    cb.TickLabels = {'0' '0.5' '1'};
    set(gca, 'FontSize', 15);
    
    xlabel('Mode (to)', 'FontSize', 18);
    xticks(1:4);
    xlim([0.5, 4.5]);
    ylabel('Mode (from)', 'FontSize', 18);
    yticks(1:4);
    ylim([0.5, 4.5]);
    
    % Add matrix value labels to each cell
    [xTxt, yTxt] = ndgrid(1:4, 1:4);
    text(xTxt(:), yTxt(:), num2cell(reshape(round(data*1e4)./1e4,1,[])'), ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
        'FontSize', 15, 'Color', 'w');
    
    xline(0.5:4.5);
    yline(0.5:4.5);
end

% Save figure
save_path = 'D:\Dropbox\Moonbrainlab\Manuscripts\Relative Phase Project\Supple\HBN';
file_name = '6_100ms_transition_unweighted';

svg_file_path = fullfile(save_path, [file_name, '.svg']);
fig_file_path = fullfile(save_path, [file_name, '.fig']);

fig = gcf;
print(fig, svg_file_path, '-dsvg');
savefig(fig, fig_file_path);

cd (save_path)
