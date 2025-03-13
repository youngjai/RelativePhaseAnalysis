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
% filter = 2; %  1- 4 (delta)
% filter = 3; %  4- 8 (theta)
% filter = 4; %  8-12 (alpha)
% bands(end+1,:) = [12, 20]; filter = length(bands); % 12-20 (low beta)
% bands(end+1,:) = [20, 30]; filter = length(bands); % 20-30 (high beta)
bands(end+1,:) = [30, 50]; filter = length(bands); % 30-50 (low gamma)
band_name = ['band_[' num2str(bands(filter,1)) '-' ...
    num2str(bands(filter,2)) ']'];
subj = subjname{idx};

% create directory to save output files
output_path = ['../preprocessing/completed_240715/' band_name '/'];
if ~exist(output_path, 'dir')
    mkdir(output_path);
end
if ~exist([output_path 'figures/'], 'dir')
    mkdir([output_path 'figures/']);
end

% eeglab run
eeglab nogui;

% eeglab load
notch = 60;
cutoff_low = 0.5; cutoff_high = 100;
load([output_path '../' subj '_wo_badchan_and_wo_reref.mat']);
n_state = length(state_list);

for st = 1:n_state
    state = state_list{st};
    data = data_state{st}';
    if ~isempty(data)
        EEG.etc.eeglabvers = '2022.1'; 
        EEG = pop_importdata('dataformat','array', 'data','data', 'srate',500);
        
        % channel location (only use 124 channels)
        EEG.chanlocs = chanlocs{st};
        EEG = eeg_checkset(EEG);
    
        EEG = pop_eegfiltnew(EEG, 'locutoff',bands(filter,1), ...
            'hicutoff',bands(filter,2));

        data_state{st} = EEG.data';
    end
end

%% save preprocessed EEG data
save_file_name = [subj '_wo_badchan_and_wo_reref_and_filtered'];
save([output_path save_file_name '.mat'], 'data_state', 'chanlocs', ...
    'state_list', '-v7.3');

%% draw power spectrogram to check the result
state_label_list = {'EO', 'EC', 'LOC', 'Anesthesia', ...
    'Burst', 'Suppression', 'Deep anesthesia', 'ROC'};
dur_list = cumsum([0, 5, 5, 3, 5, 5, 5, 5, 5]*60);
state_loc = ~cellfun(@isempty,data_state);

time_resol = 2; % Time resolution (s) = Time window size. Set as 2 s
L = 50;         % Overlap between windows (%)
lk = 0.85;      % Leakage factor for a Kaiser window : lk=0.85 is equilvalent to Appx. Hann window / Hamming window


fig = figure(1);    % make frame without display
fig.Position = [50 50 1200 300];   
clf;
%     sgtitle([subjname{idx} '_' labels{idx} '_ch_E' num2str(ch) '(Pz)'], ...
sgtitle([subjname{idx} ' ' labels{idx} ' ch_E62'], ...
    'FontSize',12, 'Interpreter','none', 'FontWeight','bold');
subplot(1,1,1);
hold on;
for st = 1:n_state
    if state_loc(st) > 0
        ch = find(ismember({chanlocs{st}.labels},'E62')); % E62 : Pz
        if isempty(ch)
            [~,ch] = min(abs([chanlocs{st}.urchan]-3-62)); % find nearest ch from Pz
            text(mean([dur_list(st) dur_list(st+1)]),-1.2, ...
                ['( ' chanlocs{st}(ch).labels ' )'], ...
                'FontSize',6, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'Color','r');
        end
        data_temp = band_pass_filter(data_state{st}(:,ch), fs, 0.5, 100); % detrending
        [p,f,t] = pspectrum(data_temp, fs, "spectrogram", ...
            TimeResolution=time_resol, OverlapPercent=L, Leakage=lk, ...
            FrequencyLimits=[1 70]);
        imagesc(t+dur_list(st),f,pow2db(p));
    end
end
set(gca,'YDir','normal');
set(gca,'FontSize',12);
xlim(dur_list([1,end]));
xticks(dur_list);
xticklabels(dur_list./60);
xlabel('Time (min.)', 'FontSize',15);
ylim([1 70]);
ylabel('Frequency (Hz)', 'FontSize',15);
cb = colorbar;
cb.Label.String = 'Power (dB)';
colormap('jet');
caxis([-25 15]); % Set colormap limits ( caxis renamed from in R2022a)
xline(dur_list, ...
    'LineWidth',1.5, 'Color','k');
text(dur_list(1:end-1)+diff(dur_list).*0.5,72*ones([1,n_state]), ...
    state_label_list, 'FontSize',6, 'FontWeight','bold', ...
    'HorizontalAlignment','center', 'Color','k');

exportgraphics(fig, [output_path 'figures/' subjname{idx} '_' labels{idx} ...
    '_pspectrum_' band_name '.png']);
