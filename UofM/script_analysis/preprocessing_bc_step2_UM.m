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

%% bad channel detection made by Younghwa CHA

% clear, close all;

load('../results_yjp/UM_info_ch128.mat');
% idx = 1;
subj = subjname{idx};

% create directory to save output files
output_path = '../preprocessing/completed_240715/';
if ~exist(output_path, 'dir')
    mkdir(output_path);
end
if ~exist([output_path 'figures/'], 'dir')
    mkdir([output_path 'figures/']);
end

% load EGI-GSN 129 channel information (E129: Cz)
load('references/ch_129.mat');
load('references/EGI_GSN_129.mat');
ch_128_labels = {EGI_GSN_128.labels};

criterion_p = 0.25;
ch_file_name = ['common_wb_124_bc_idx_crit' ...
    'erion_' num2str(criterion_p)];
ch_file_name = strrep(ch_file_name, '.', 'p');
load([output_path '../badchan_240715/' ch_file_name '.mat'], 'ch_loc_common');

%% make an index for removing bad channels
good_session_thr = 124*0.75; % bad session : 25 percent of 124 channels

% eeglab run
eeglab nogui;

% eeglab load
notch = 60;
cutoff_low = 0.5; cutoff_high = 100;
load_file_name = [subj '_notch_' num2str(notch) ...
    '_bandpass_' num2str(cutoff_low) '_' num2str(cutoff_high)];
load_file_name = strrep(load_file_name, '.', 'p');
load([output_path '../notch_bandpass_240715/' load_file_name '.mat']);
n_state = length(state_list);
chanlocs = cell([1,n_state]);

for st = 1:n_state
    state = state_list{st};
    data = data_state{st}';
    if ~isempty(data)
        EEG.etc.eeglabvers = '2022.1'; 
        EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
        
        % channel location (only use 124 channels)
        EEG.chanlocs = EGI_GSN_128(1:124);
        EEG = eeg_checkset(EEG);
        
        % remove common bad channels
        ch_labels = {ch_loc_common.labels};
        EEG = pop_select(EEG, 'channel', ch_labels);
        EEG = eeg_checkset(EEG);
        
        % re-remove bad channels
        trim_low = 1e-4;    trim_high = 100;
        EEG1 = trimOutlier(EEG, trim_low, trim_high, Inf, 0);
        EEG1 = eeg_checkset(EEG1);
        
        bad_ch_info(st).good = EEG1.chanlocs;
        try 
          bad_ch_info(st).bad = EEG1.chaninfo.removedchans;
        catch
          warning('No channel removed.');
          bad_ch_info(st).bad = [];
        end
        bad_ch_info(st).name = subj;
        clear EEG1;
            
        %% Make a common bad channel list
        ch_labels_idx = [bad_ch_info(st).bad.urchan]-3; % due to urchan starts from 4
        
        % find symmetric channels of bc
        % before loading, plz make a symmetric array considering your eeg channal map
        
        bad_sym = ch_129(2, ch_labels_idx);
        bad_combined = unique([bad_sym ch_labels_idx]);
        idx_remove = find(bad_combined == 999);
        bad_combined(idx_remove) = [];
        
        bc_sym = bad_combined;
        good_channel = setdiff(1:124, bc_sym);
        
        ch_labels = ch_128_labels(good_channel);

        if length(good_channel) >= good_session_thr
            chanlocs{st} = EGI_GSN_128(good_channel);
            
            % remove common bad channels
            EEG = pop_select(EEG, 'channel', ch_labels);
            EEG = eeg_checkset(EEG);
            
%             % rereferencing using average reference (that's what [] means)
%             EEG = pop_reref(EEG, []);
        
            data_state{st} = EEG.data';
        else
            data_state{st} = [];
        end
        clear EEG;

        fig = figure(1);
        fig.Position = [100 100 1500 450];
        clf;
        sgtitle(sprintf('%s, state: %s\nbad:%d/129, good:%d/129', ...
            subjname{idx}, state, 129-length(good_channel), length(good_channel)), ...
            'Interpreter','none', 'FontWeight','bold');
        subplot(1,3,1);
        hold on; grid on;
        plot3(-[EGI_GSN_128([125:129 bc_sym]).Y], ...
            [EGI_GSN_128([125:129 bc_sym]).X], ...
            [EGI_GSN_128([125:129 bc_sym]).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
        plot3(-[EGI_GSN_128(good_channel).Y], ...
            [EGI_GSN_128(good_channel).X], ...
            [EGI_GSN_128(good_channel).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
        text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
            {EGI_GSN_128.labels}, 'FontSize',6, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');        
        view(0,90);
        
        subplot(1,3,2);
        hold on; grid on;
        plot3(-[EGI_GSN_128([125:129 bc_sym]).Y], ...
            [EGI_GSN_128([125:129 bc_sym]).X], ...
            [EGI_GSN_128([125:129 bc_sym]).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
        plot3(-[EGI_GSN_128(good_channel).Y], ...
            [EGI_GSN_128(good_channel).X], ...
            [EGI_GSN_128(good_channel).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
        text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
            {EGI_GSN_128.labels}, 'FontSize',6, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');        
        view(-91,2);
        
        subplot(1,3,3);
        hold on; grid on;
        plot3(-[EGI_GSN_128([125:129 bc_sym]).Y], ...
            [EGI_GSN_128([125:129 bc_sym]).X], ...
            [EGI_GSN_128([125:129 bc_sym]).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#EDB120');
        plot3(-[EGI_GSN_128(good_channel).Y], ...
            [EGI_GSN_128(good_channel).X], ...
            [EGI_GSN_128(good_channel).Z], 'o', 'MarkerSize',4, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','#0072BD');
        text(-[EGI_GSN_128.Y], [EGI_GSN_128.X], [EGI_GSN_128.Z], ...
            {EGI_GSN_128.labels}, 'FontSize',6, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');        
        view(-180,-1);
        
        exportgraphics(fig, [output_path 'figures/' subjname{idx} ...
            '_st_' state '_wo_reref_badchan.png']);
    end
end

%% save preprocessed EEG data
save([output_path 'bad_ch_info_' subj '_re.mat'], 'bad_ch_info');

save_file_name = [subj '_wo_badchan_and_wo_reref'];
save([output_path save_file_name '.mat'], 'data_state', 'chanlocs', ...
    'cutoff_low', 'cutoff_high', 'state_list', '-v7.3');
