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
% idx = 6;
subj = subjname{idx};

% create directory to save output files
output_path = '../preprocessing/badchan_240715/';
if ~exist(output_path, 'dir')
    mkdir(output_path);
end
if ~exist([output_path 'figures/'], 'dir')
    mkdir([output_path 'figures/']);
end

% load EGI-GSN 129 channel information (E129: Cz)
load('references/EGI_GSN_129.mat');

%% load notchfiltered EEG data
notch = 60;
cutoff_low = 0.5; cutoff_high = 100;
save_file_name = [subj '_notch_' num2str(notch) ...
    '_bandpass_' num2str(cutoff_low) '_' num2str(cutoff_high)];
save_file_name = strrep(save_file_name, '.', 'p');
load(['../preprocessing/notch_bandpass_240715/' save_file_name '.mat']);
n_state = length(state_list);

% eeglab load
eeglab nogui;

% in matdata500, 'data', 'chanlocs', 'etime'
% data = data_state{st}';
% EEG.etc.eeglabvers = '2022.1'; 
% % this tracks which version of EEGLAB is being used, you may ignore it
% EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
% EEG.setname = subj;
% EEG = eeg_checkset(EEG);

%% find bad channel (not save EEG at this time)
% trim_low = 1e-4;    trim_high = 100;
% % function EEG = trimOutlier(EEG, channelSdLowerBound, ...
% % channelSdUpperBound, amplitudeThreshold, pointSpreadWidth)
% EEG = trimOutlier(EEG, trim_low, trim_high, Inf, 0);
% EEG = eeg_checkset(EEG);
% 
% bad_ch_info.good = EEG.chanlocs;
% try 
%   bad_ch_info.bad = EEG.chaninfo.removedchans;
% catch
%   warning('No channel removed.');
%   bad_ch_info.bad = [];
% end
% 
% bad_ch_info.name = subj;
% 
% save([output_path 'bad_ch_info_' subj '.mat'], 'bad_ch_info');

%% find bad channel (not save EEG at this time) tmp
trim_low = 1e-4;    trim_high = 100;
% function EEG = trimOutlier(EEG, channelSdLowerBound, ...
% channelSdUpperBound, amplitudeThreshold, pointSpreadWidth)

for st = 1:n_state
    state = state_list{st};
    data = data_state{st}';
    if ~isempty(data)
        EEG.etc.eeglabvers = '2022.1'; 
        % this tracks which version of EEGLAB is being used, you may ignore it
        EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
        EEG.setname = subj;
        EEG = eeg_checkset(EEG);

        % channel location (only use 124 channels)
        EEG.chanlocs = EGI_GSN_128(1:124);
        ch_labels = {EGI_GSN_128(1:124).labels};
        
        % remove common bad channels
        EEG = pop_select(EEG, 'channel', ch_labels);
        EEG = eeg_checkset(EEG);

        %% find bad channel (not save EEG at this time)
        EEG = trimOutlier(EEG, trim_low, trim_high, Inf, 0);
        EEG = eeg_checkset(EEG);
        
        bad_ch_info(st).good = EEG.chanlocs;
        try 
          bad_ch_info(st).bad = EEG.chaninfo.removedchans;
        catch
          warning('No channel removed.');
          bad_ch_info(st).bad = [];
        end
    
        bad_ch_info(st).name = [subj ' | trim_low: ' num2str(trim_low) ...
            ' trim_high: ' num2str(trim_high)];
        
    end
    save([output_path 'bad_ch_info_' subj '.mat'], 'bad_ch_info');
end

%% find bad channel (not save EEG at this time) tmp
% trim_low = 1e-4;    trim_high = 100;
% % function EEG = trimOutlier(EEG, channelSdLowerBound, ...
% % channelSdUpperBound, amplitudeThreshold, pointSpreadWidth)
% 
% trim_high_list = [100 150 200];
% 
% n_state = length(state_list);
% for st = 1:n_state
%     state = state_list{st};
%     data = data_state{st}';
%     if ~isempty(data)
%         EEG.etc.eeglabvers = '2022.1'; 
%         % this tracks which version of EEGLAB is being used, you may ignore it
%         EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
%         EEG.setname = subj;
%         EEG = eeg_checkset(EEG);
% 
%         % channel location (only use 124 channels)
%         EEG.chanlocs = EGI_GSN_128(1:124);
%         ch_labels = {EGI_GSN_128(1:124).labels};
%         
%         % remove common bad channels
%         EEG = pop_select(EEG, 'channel', ch_labels);
%         EEG = eeg_checkset(EEG);
% 
%         for tr = 1:length(trim_high_list)
%             trim_high = trim_high_list(tr);
%             EEG1 = trimOutlier(EEG, trim_low, trim_high, Inf, 0);
%             EEG1 = eeg_checkset(EEG1);
%             
%             bad_ch_info(tr).good = EEG1.chanlocs;
%             try 
%               bad_ch_info(tr).bad = EEG1.chaninfo.removedchans;
%             catch
%               warning('No channel removed.');
%               bad_ch_info(tr).bad = [];
%             end
%         
%             bad_ch_info(tr).name = [subj ' | trim_low: ' num2str(trim_low) ...
%                 ' trim_high: ' num2str(trim_high)];
%         end
%         
%         save([output_path 'bad_ch_info_' subj '_st_' state '.mat'], 'bad_ch_info');
%     end
% end
