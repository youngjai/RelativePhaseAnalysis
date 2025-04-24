%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load calculated relative phase results
% 2. Compute group level averages for the results
% 3. Merge features from Control and ADHD groups for comparison
% 4. Load behavioral scores (e.g., SWAN) and diagnostic labels
% 5. Build logistic regression models for classification 
% 6. Correlate predicted probabilities with SWAN scores

% 2025. 03. 27. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc;

%% data extract

cd('D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result')
load("C_ec_cal_sub")
load("C_eo_cal_sub")
load("Ain_ec_cal_sub")
load("Ain_eo_cal_sub")

condition = 'alpha';  % 'theta' or 'delta'

files = {
    ['C_ec'],
    ['C_eo'],
    ['Ain_ec'],
    ['Ain_eo']
};

% mode and dwell 
for i = 1:length(files)
    file_name = files{i};
    file_task = file_name;

    var_name = ['prop_sub_' file_task ];       


    eval([' prop_sub =' var_name ';']);  
    
    for subj_no = 1:length(prop_sub)
        
        data = prop_sub(subj_no, 1:5);
    
        for j = 1:5
        
            mode(j, 1:4) = [data(j).occurrence];
            dwell(j, 1:4) = [data(j).dwell_time];

            for z = 1:4
                trans(j, z) = [data(j).trans_freq_mat(z, 4)];
                trans_r(j, z) = [data(j).trans_freq_mat(4, z)];
            end
    
        end
    
        sub_mode(subj_no, 1:4) = mean(mode);
        sub_dwell(subj_no, 1:4) = mean(dwell);
        sub_trans(subj_no, :) = mean(trans);
%         sub_trans(subj_no, :) = mean(trans_r);
    
    end

    var_name_mode = ['sub_mode_' file_task]; 
    var_name_dwell = ['sub_dwell_' file_task]; 
    var_name_trans = ['sub_trans_' file_task]; 

    eval([var_name_mode ' = sub_mode;']);  
    eval([var_name_dwell ' = sub_dwell;']);  
    eval([var_name_trans ' = sub_trans;']);  

    mean_mode = mean(sub_mode); % group
    mean_dwell = mean(sub_dwell); % group
    mean_trans = mean(sub_trans); % group

    se_mode = std(sub_mode)/sqrt(length(sub_mode));
    se_dwell = std(sub_dwell)/sqrt(length(sub_dwell));
    se_trans = std(sub_trans)/sqrt(length(sub_trans));

    var_name_mode_m = ['mean_mode_' file_task]; 
    var_name_dwell_m = ['mean_dwell_' file_task]; 
    var_name_trans_m = ['mean_trans_' file_task]; 

    eval([var_name_mode_m ' = mean_mode;']);  
    eval([var_name_dwell_m ' = mean_dwell;']);  
    eval([var_name_trans_m ' = mean_trans;']);  

    var_name_mode_se = ['se_mode_' file_task]; 
    var_name_dwell_se = ['se_dwell_' file_task]; 
    var_name_trans_se = ['se_trans_' file_task]; 

    eval([var_name_mode_se ' = se_mode;']);  
    eval([var_name_dwell_se ' = se_dwell;']);  
    eval([var_name_trans_se ' = se_trans;']);  

    clear sub_mode sub_dwell sub_trans mode dwell trans trans_r data w_prop_sub mean_mode mean_dwell mean_trans se_mode se_dwell se_trans

end


%% variables

total = readtable('C:\Users\user\Downloads\HBN_Rating_total_sorting.csv');

% labels for diagnosing
labels(1:66) = 0;
labels(67:106) = 1;
labels = labels';

% swan score
swan = total.SWAN_IN_Avg_3;

ec_m_b(1:66, :) = sub_mode_C_ec;
ec_m_b(67:106, :) = sub_mode_Ain_ec;

eo_m_b(1:66, :) = sub_mode_C_eo;
eo_m_b(67:106, :) = sub_mode_Ain_eo;

ec_d_b(1:66, :) = sub_dwell_C_ec;
ec_d_b(67:106, :) = sub_dwell_Ain_ec;

eo_d_b(1:66, :) = sub_dwell_C_eo;
eo_d_b(67:106, :) = sub_dwell_Ain_eo;

ec_t_b(1:66, :) = sub_trans_C_ec;
ec_t_b(67:106, :) = sub_trans_Ain_ec;

eo_t_b(1:66, :) = sub_trans_C_eo;
eo_t_b(67:106, :) = sub_trans_Ain_eo;


% final combination (t-test star)
X = [ec_m_b(:, 4), eo_m_b(:, 4), ec_d_b(:, 4), eo_d_b(:, 4), ec_t_b(:, 4), eo_t_b(:, 4)];
X_zscored = zscore(X);

% tbr combination 
X_tbr = [fz_eo, fz_ec, pz_eo, pz_ec];
X_tbr_zscored = zscore(X_tbr);

y = labels;

[beta, dev, stats] = glmfit(X, y, 'binomial', 'link', 'logit');
[beta_tbr, dev_tbr, stats_tbr] = glmfit(X_tbr, y, 'binomial', 'link', 'logit');

% calculation for prediction values based on logistic model
p_hat= glmval(beta, X, 'logit');
p_hat_tbr= glmval(beta_tbr, X_tbr, 'logit');

% ROC
[X_ROC, Y_ROC, T_ROC, AUC_ROC] = perfcurve(y, p_hat, 1);
[X_ROC_tbr, Y_ROC_tbr, T_ROC_tbr, AUC_ROC_tbr] = perfcurve(y, p_hat_tbr, 1);

% % for checking
% figure;
% plot(X_ROC, Y_ROC);
% xlabel('False Positive Rate');
% ylabel('True Positive Rate');
% title('ROC Curve');
% grid on;


%% corr survey
[coef1, pval1] = corr(p_hat, swan, 'Rows', 'complete');

fig = figure(101);
hold on;

% first group (subjects 1 to 66) - green dots
scatter(p_hat(1:66), swan(1:66), 'MarkerEdgeColor', [27/255, 94/255, 32/255], ...
    'MarkerFaceColor', [27/255, 94/255, 32/255]);

% second group (subjects 67 to 106) - red dots
scatter(p_hat(67:106), swan(67:106), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% add regression line (correlation line)
% fit a linear regression (1st-degree polynomial)
coeffs = polyfit(p_hat, swan, 1);

x_fit = linspace(min(p_hat), max(p_hat), 100);

% calculate corresponding y values from the fitted line
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);

xlabel('RPA Logistic Model', 'FontSize', 14);
ylabel('SWAN-Inattentive', 'FontSize', 14);
title('Correlation between RPA and SWAN', 'FontSize', 16);

legend({'Control', 'ADHD'}, 'FontSize', 12, 'Location', 'southeast');

grid on;
hold off;


