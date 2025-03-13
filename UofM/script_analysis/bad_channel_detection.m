%% bad channels detection
% % made by youngjai 2023.05.02.
% 
% load_path = [pwd() '/..']; % Project location
% save_path = [load_path '/preprocessing'];
% 
% % 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
% load([load_path '/results_yjp/UM_info_AE.mat']);
% % 'UM_*.mat' contain 'data', 'chanlocs', 'etime' (equivalent to eventtime), and 'labels'
% 
% outliers = zeros(n_sub_ae,n_ch);
% for idx = 1:n_sub_ae
%     % load([load_path '/matdata500/' subjname_ae{idx} '.mat']);
%     % LBF = 1; HBF = 50;
%     % data_bp = band_pass_filter(data, fs, LBF, HBF);
%     % clear data;
%     
%     % size_ch     = size(chanlocs,2); % the number of channels
%     % chan_idx    = 1:size_ch;
%     % [data_bp_rf] = reref_data(double(data_bp), 'all', chan_idx);
%     % clear data_bp;
% 
%     % k = kurtosis(data_bp_rf);
% 
%     load([save_path '/slicing/band_bf/' subjname_ae{idx} '_all_' labels_ae{idx} '_band_bf_prepr.mat']);
%     k = kurtosis(data);
% 
%     outliers(idx,:) = isoutlier(k);
% end
% bad_chan = sum(outliers);
% 
% % save([save_path '/UM_bad_channels_info.mat'], 'fs', 'n_sub_ae', 'labels_ae', 'n_ch', 'bad_chan', 'outliers');
% save([save_path '/UM_bad_channels_info1.mat'], 'fs', 'n_sub_ae', 'labels_ae', 'n_ch', 'bad_chan', 'outliers');

%%
load_path = [pwd() '/..']; % Project location
save_path = [load_path '/preprocessing'];
% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
load([load_path '/results_yjp/UM_info_AE.mat']);
% load([save_path '/UM_bad_channels_info.mat']);
load([save_path '/UM_bad_channels_info1.mat']);

fig = figure('NumberTitle','off', 'Visible','off');    % make frame without display
fig.Position(3:4) = [1200,600]; 
% figure(23);
clf; 
% idx = 1;

% load([load_path '/matdata500/' subjname_ae{idx} '.mat']);
% 
% LBF = 1; HBF = 50;
% data_bp = band_pass_filter(data, fs, LBF, HBF);
% clear data;
% 
% size_ch     = size(chanlocs,2); % the number of channels
% chan_idx    = 1:size_ch;
% [data_bp_rf] = reref_data(double(data_bp), 'all', chan_idx);
% clear data_bp;
% data = data_bp_rf;
% clear data_bp_rf;

load([save_path '/slicing/band_bf/' subjname_ae{idx} '_all_' labels_ae{idx} '_band_bf_prepr.mat']);

bad_chan = logical(outliers);

subplot(1,3,1);
% Caution! to change the X and Y axis!!
X = [chanlocs.Y];
Y = [chanlocs.X];
Z = [chanlocs.Z];
% plot3(X,Y,Z,'ko','MarkerFaceColor','b');
plot(X,Y,'ko','MarkerFaceColor','b');
hold on;
grid on;
% plot3(X(bad_chan(idx,:)),Y(bad_chan(idx,:)),Z(bad_chan(idx,:)),'ko','MarkerFaceColor','g');
plot(X(bad_chan(idx,:)),Y(bad_chan(idx,:)),'ko','MarkerFaceColor','g');

if sum(bad_chan(idx,:))
    subplot(1,3,2);
    pspectrum(data(:,bad_chan(idx,:)),fs,'FrequencyLimits',[1 50]);
end

subplot(1,3,3);
pspectrum(data(:,~bad_chan(idx,:)),fs,'FrequencyLimits',[1 50]);

% exportgraphics(fig, [load_path '/results_yjp/bad_channel/' subjname_ae{idx} '_bad_channel.png'], ...
%     'Resolution',600)
exportgraphics(fig, [load_path '/results_yjp/bad_channel1/' subjname_ae{idx} '_bad_channel.png'], ...
    'Resolution',600)
