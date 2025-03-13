%% Preprocessing Pipeline
% This code is based on Younghwa's preprocessing code
% in \Dropbox\Moonbrainlab\HBN\code\HBN_wb_preprocessing_resting_new.m
%
% Step 1. Do slicing
%
% Step 2. Do Notch filtering
% HBN, UofM notch filter: 60 Hz
%
% Step 3. Do Band-pass filtering
% Bandpass filter: 0.5-100 Hz
%
% Step 4. Find a comman bad channel list
% Find a common bad channel list by checking SD of channel amplitude rejection threshold (10-4μV, 100μV)
% a. HBN: number of bad channel criteria, 25% of all subjects.
%
% Step 5. Remove bad channels
% Remove bad channels according to SD without the common bad channel list
% a. If the number of bad channels with over 25% of the number of base channels, the subject will be excluded 
%
% Step 6. Re-reference using average referencing
%
% Step 7. Bandpass filter for each band
% 1-4(delta), 4-8(theta), 8-12(alpha), 13-30(beta), 30-50(low-gamma)


% bad channel detection made by Younghwa CHA

%% Make a common bad channel list
clear, close all;

load('../results_yjp/UM_info_ch128.mat');
output_path = '../preprocessing/badchan_240715/';

% state_list = {'EO', 'EC', 'AI1', 'LOC', 'AI2', 'AN', 'BS', 'DS', 'DA', 'ROC'};
state_list = {'EO', 'EC', 'LOC', 'AN', 'BS', 'DS', 'DA', 'ROC'};
n_state = length(state_list);

% change bad index name and load
for idx = 1:n_sub
    load([output_path 'bad_ch_info_' subjname{idx} '.mat']);

    for st = 1:size(bad_ch_info,2)
        bc_idx(idx,st).good = bad_ch_info(st).good;
        bc_idx(idx,st).bad = bad_ch_info(st).bad;
        bc_idx(idx,st).name = bad_ch_info(st).name;
    end
end

% make a common bad channel list
n_session = 0;
for idx = 1:n_sub
    for st = 1:n_state
        if ~isempty(bc_idx(idx,st).bad)
            urchan = [bc_idx(idx,st).bad.urchan];
            n_session = n_session + 1;
            bc_list{n_session} = urchan;
            clear('urchan');
        end 
    end
end


% make a count array for bad channels
count_array = zeros(1, 132);

% update count array
for idx = 1:n_session
    temp = bc_list{idx};
    if ~isempty(temp)
        for i = 1:length(temp)
            ch = temp(i);
            count_array(ch) = count_array(ch) + 1;
        end
    end
end

% criterion: 25% of the total number of the data
criterion_p = 0.25;
criterion = n_session*criterion_p;

% fit to 129 (including cz)
count_array([1 2 3]) = [];

% common bad channels
common_bc = find(count_array > criterion);

% find symmetric channels of common bc
% before loading, plz make a symmetric array considering your eeg channal map
load('references/ch_129.mat');

bad_sym = ch_129(2, common_bc);
bad_combined = unique([bad_sym common_bc]);
idx_remove = find(bad_combined == 999);
bad_combined(idx_remove) = [];

common_bc_sym = bad_combined;
good_channel = setdiff(1:124, common_bc_sym);

save_file_name = ['common_wb_124_bc_idx_criterion_' num2str(criterion_p)];
save_file_name = strrep(save_file_name, '.', 'p');
save([output_path save_file_name '.mat'], 'common_*', 'good_channel', '-v7.3');


%% Remove common bad channels
% load the common bad channel list
load([output_path save_file_name '.mat']);

% make a new channel location
% load previously saved channel information file
load('references/EGI_GSN_129.mat');

EGI_GSN_128 = EGI_GSN_128(good_channel);

% channel location without common bad channels
ch_loc_common = EGI_GSN_128;

% save([output_path 'common_ch_loc_info.mat'], 'ch_loc_common');
save([output_path save_file_name '.mat'], 'ch_loc_common', '-append');

%% draw comman bad channels
load('references/EGI_GSN_129.mat');

fig = figure(1);
fig.Position = [100 100 1500 450];
clf;
sgtitle(sprintf('UofM general anesthesia 128ch location\nbad:%d/129, good:%d/129', ...
    129-length(ch_loc_common), length(ch_loc_common)), ...
    'FontWeight','bold');
subplot(1,3,1);
hold on; grid on;
plot3(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
    {EGI_GSN_128.labels}, 'FontSize',6, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');

plot3(-[ch_loc_common.Y], [ch_loc_common.X], [ch_loc_common.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
view(0,90);

subplot(1,3,2);
hold on; grid on;
plot3(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
    {EGI_GSN_128.labels}, 'FontSize',6, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');

plot3(-[ch_loc_common.Y], [ch_loc_common.X], [ch_loc_common.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
view(-91,2);

subplot(1,3,3);
hold on; grid on;
plot3(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
    {EGI_GSN_128.labels}, 'FontSize',6, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');

plot3(-[ch_loc_common.Y], [ch_loc_common.X], [ch_loc_common.Z], 'o', 'MarkerSize',4, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
view(-180,-1);

exportgraphics(fig, [output_path '/figures/' save_file_name '.png']);
