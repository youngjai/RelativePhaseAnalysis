%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load calculated relative phase results
% 2. Compute group-level averages from the results
% 3. Merge features from Control and ADHD groups for comparison
% 4. Load behavioral scores (e.g., SWAN) and diagnostic labels
% 5. Build logistic regression models for classification 
% 6. Calculate ROC curves 
% 7. Perform bootstrapping
% 8. Plot figures

% 2024.10.14 - Younghwa Cha
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


%% ROC
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


%% bootstraping
% Generate example data (replace with actual data)
n = 106;
p = size(X, 2);

% Number of bootstrap iterations
num_bootstrap = 1000;

% Initialize array to store AUC values
auc_scores = zeros(num_bootstrap, 1);
auc_scores_tbr = zeros(num_bootstrap, 1);

for i = 1:num_bootstrap
    % Bootstrap sampling (sampling with replacement)
    idx = randsample(n, n, true);
    X_bootstrap = X(idx, :);
    X_tbr_bootstrap = X_tbr(idx, :);
    y_bootstrap = y(idx);
    
    % Z-score normalization
    % X_bootstrap_z = zscore(X_bootstrap);
    
    % Train logistic regression model
    [beta_bootstrap, ~, ~] = glmfit(X_bootstrap, y_bootstrap, 'binomial', 'link', 'logit');
    [beta_tbr_bootstrap, ~, ~] = glmfit(X_tbr_bootstrap, y_bootstrap, 'binomial', 'link', 'logit');
    
    % Predict probabilities on the original dataset
    X_z = (X - mean(X_bootstrap)) ./ std(X_bootstrap); % Normalize original data
    p_hat_bootstrap = glmval(beta_bootstrap, X_z, 'logit');
    X_tbr_z = (X_tbr - mean(X_tbr_bootstrap)) ./ std(X_tbr_bootstrap); % Normalize original data
    p_hat_tbr_bootstrap = glmval(beta_tbr_bootstrap, X_tbr_z, 'logit');
    
    % Compute ROC curve
    [X_b, Y_b, ~, auc_b] = perfcurve(y, p_hat_bootstrap, 1);
    auc_scores(i) = auc_b;
    [X_b_tbr, Y_b_tbr, ~, auc_b_tbr] = perfcurve(y, p_hat_tbr_bootstrap, 1);
    auc_scores_tbr(i) = auc_b_tbr;

    X_b_all{i} = X_b;
    Y_b_all{i} = Y_b;

    X_b_tbr_all{i} = X_b_tbr;
    Y_b_tbr_all{i} = Y_b_tbr;
end

% Calculate mean and standard deviation of AUC
mean_auc = mean(auc_scores);
std_auc = std(auc_scores);
mean_auc_tbr = mean(auc_scores_tbr);
std_auc_tbr = std(auc_scores_tbr);

% Initialize for averaged ROC curve (removing N/A)
mean_X = zeros(107, 1); 
mean_Y = zeros(107, 1);

mean_X_tbr = zeros(107, 1); 
mean_Y_tbr = zeros(107, 1);

for i = 1:length(X_b_all)
    mean_X = mean_X + X_b_all{i};  % Accumulate data from each cell
    mean_Y = mean_Y + Y_b_all{i};
    mean_X_tbr = mean_X_tbr + X_b_tbr_all{i};
    mean_Y_tbr = mean_Y_tbr + Y_b_tbr_all{i};
end

% Final averaged ROC curve
mean_X = mean_X / length(X_b_all); 
mean_Y = mean_Y / length(Y_b_all); 
mean_X_tbr = mean_X_tbr / length(X_b_tbr_all); 
mean_Y_tbr = mean_Y_tbr / length(Y_b_tbr_all); 


%% figure
colors(1, :) = [0.8500 0.3250 0.0980];
colors(2, :) = [0.9290 0.6940 0.1250];
colors(3, :) = [0 0.4470 0.7410];
colors(4, :) = [0.3010 0.7450 0.9330];

figure(1)
hold on
grid on
% RP (X_ROC, Y_ROC)
plot(X_ROC, Y_ROC, 'Color', colors(1, :), 'LineWidth', 1.5);

% RP_Boots (mean_X_s, mean_Y_s)
plot(mean_X, mean_Y, 'Color', colors(2, :), 'LineWidth', 1.5);

% TBR (X_ROC_t, Y_ROC_t)
plot(X_ROC_tbr, Y_ROC_tbr, 'Color', colors(3, :), 'LineWidth', 1.5);

% TBR_Boots (mean_X_t, mean_Y_t)
plot(mean_X_tbr, mean_Y_tbr, 'Color', colors(4, :), 'LineWidth', 1.5);


plot([0, 1], [0, 1], '--', 'Color', [0.5, 0.5, 0.5]); 

xlabel('False Positive Rate', 'FontSize', 24);
ylabel('True Positive Rate', 'FontSize', 24);
title('ROC Curve Combination', 'FontSize', 24);
grid on;

legend('RP','RP Boots', 'TBR', 'TBR Boots','criteria', 'location', 'southeast')


