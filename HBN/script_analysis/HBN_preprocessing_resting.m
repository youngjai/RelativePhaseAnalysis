%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load EEG files
% 2. Notch filtering - 60Hz
% 3. Bandpass filtering - 0.5Hz-100Hz
% 4. Find a common bad channel list by checking SD (10-4μV, 100μV)
% 5. Remove common bad channels (all subject) 
% 6. Re-remove bad channels using SD (ch by ch)

% 2024. 01. 16. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% Load EEG files
% change file path and load subj_list
% create your data list and load it.

list = readtable('C:\Users\user\Downloads\HBN_Control.csv');
control = table2array(list(:, :));
subj_list = control;

list = readtable('C:\Users\user\Downloads\HBN_ADHDin.csv');
control = table2array(list(:, :));
subj_list2 = control;

com_subject_list = vertcat(subj_list, subj_list2);


% load previously saved channel information file
load(['D:\HealthyBrainNetwork_above11\channel_location\EGI_GSN_129.mat']);
load(['D:\HealthyBrainNetwork_above11\channel_location\EGI_GSN_124_correct.mat']);

% make a folder 
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = 'D:\HealthyBrainNetwork_above11\EEGDownload4\C_w_124_bp_original_new\';
outputDir2 = 'D:\HealthyBrainNetwork_above11\EEGDownload4\Ain_w_124_bp_original_new\';
mkdir([outputDir1])
mkdir([outputDir2])

%% Notch filter / Bandpass 
% select subj_no considering subj_list

for subj_no = 1:length(com_subject_list)
    
    % load the file
    input_fileName = char(com_subject_list(subj_no));
%     fileName = erase(input_fileName, '.mat')
    fileName = input_fileName;

    % file path
    unzip = [outputDir 'original\unzip2\'];
    file_path = [unzip com_subject_list(subj_no) '\EEG\raw\mat_format\'];
    file_path = strjoin(file_path);
    file_path = file_path(find(~isspace(file_path)));

    % state
    state = 'RestingState.mat';
    stateName = 'Resting'

    % eeglab run
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % eeglab load
    EEG = pop_loadset('filename',state,'filepath',file_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    %eeglab redraw

    % eeglab notchfilter and bandpass filtering(0.5~100hz)
    EEG = eeglab_notchFilter(EEG,outputDir,60)
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',100,'plotfreqz',1);
    close gcf

    event = EEG.event;
    event_sample = zeros(length(event), 1);
    event_type = transpose({event(:).type});
    
    for t = 1:length(event)
    event_sample(t, 1) = getfield(event,{t},'sample');
    end
    
    % erase break period(90)
    event_sample = event_sample(3:13);
    event_type = event_type(3:13);
    
    event_sample(1, 2) = 1;
    for i = 2:length(event_sample)
    event_sample(i, 2) = event_sample(i, 1)-event_sample((i-1), 1)+event_sample((i-1), 2);
    end
    
    event_time = event_sample*2/1000/60;
    
    start_point = event_sample(1, 1);
    last_point = event_sample(end, 1);
    
    
    EEG.data = EEG.data(:, start_point:last_point);
    EEG.event = EEG.event(3:13);
    data = EEG.data;
    Fs = EEG.srate;
    size_time = length(data);
    time = 0:1000/Fs:(size_time-1)*1000/Fs;
    size_ch = EEG.nbchan-5;
    EEG.times = time;
    EEG.pnts = length(data);
    
    % change channel number
    EEG.data = EEG.data(1:124, :);
    EEG.nbchan = 124;
    EEG.nbchans = 124;

    % channel location
    EEG.chanlocs = EGI_GSN_124;

    if subj_no <= length(subj_list)
        save_path_bp = [outputDir1 fileName '_bp'] %bp means bandpass filtering
    else 
        save_path_bp = [outputDir2 fileName '_bp'] 
    end
    
    save(save_path_bp , '-v7.3', 'EEG', 'fileName')


    % find bad channel (not save EEG at this time)

    EEG = trimOutlier(EEG, 1e-4, 100, Inf, 0);
    EEG = eeg_checkset( EEG );

    bc_idx(subj_no).good = EEG.chanlocs;

    try 
      bc_idx(subj_no).bad = EEG.chaninfo.removedchans;
    catch
      warning('No channel removed.');
      bc_idx(subj_no).bad = [];
    end

    bc_idx(subj_no).name = fileName

    % DO NOT SAVE THE EEG FILE

end

save([outputDir 'bc_idx_wb_124_C_Ain_original'], 'bc_idx')


%% Make a common bad channel list

clear; clc;

% change file path 
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';

% change bad ch index name and load
load([outputDir 'bc_idx_wb_124_C_Ain_original'])

% make a common bad channel list
for subj_no = 1:length(bc_idx)

    if length(bc_idx(subj_no).bad) ~= 0 
        for i = 1:length(bc_idx(subj_no).bad)
            urchan(i) = bc_idx(subj_no).bad(i).urchan;
        end
    
        bc_list{subj_no} = urchan;
        clear('urchan')
    end 
end

% make a count array for bad channels
count_array = zeros(1, 124);

% update count array
for subj_no = 1:length(bc_list)
    temp = bc_list{subj_no};
    if ~isempty(temp)
        for i = 1:length(temp)
            ch = temp(i);
            count_array(1, ch) = count_array(1, ch) + 1;
        end
    end
end

% criterion: 25% of the total number of the data
criterion = length(bc_list)*0.25;

% % fit to 129 (including cz)
% count_array(:, [1 2 3]) = [];

% common bad channels
common_bc = find(count_array > criterion);

% find symmetric channels of common bc
% before loading, plz make a symmetric array considering your eeg channal map
load('D:\HealthyBrainNetwork_above11\channel_location\ch_129.mat')

bad_sym = ch_129(2, common_bc);
bad_combined = unique([bad_sym common_bc]);
idx_remove = find(bad_combined == 999);
bad_combined(idx_remove) = [];

common_bc_sym = bad_combined;
good_channel = setdiff ( 1:124 , common_bc_sym );

save([outputDir 'common_wb_124_bc_idx_25p_original'], 'common_*', 'good_channel')


%% Remove common bad channels
% From this point, it should be done by group.

clear; clc;

% change file path
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = [outputDir 'Ain_w_124_bp_original_new\']
outputDir2 = [outputDir 'Ain_w_124_bp_bc_original_new\']
outputDir3 = 'D:\HealthyBrainNetwork_above11\channel_location\';
mkdir([outputDir2])

% load the common bad channel list
load(['D:\HealthyBrainNetwork_above11\EEGDownload4\common_wb_124_bc_idx_25p_original.mat']);
load('D:\HealthyBrainNetwork_above11\channel_location\EGI_GSN_124_correct.mat');

% load subj_list
subj_list = dir([outputDir1 '*_*']);

% make an index for removing bad channels (only index saving) 
for subj_no = 1:length(subj_list)
    
    cd(outputDir1)
    
    % load the file
    input_fileName = char(subj_list(subj_no).name);
    fileName = erase(input_fileName, '.mat')

    % eeglab run
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % eeglab load
    EEG = pop_loadset('filename',fileName,'filepath',outputDir1);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    %eeglab redraw

    data = EEG.data;
    Fs = EEG.srate;

    bad_channel = common_bc_sym;
    good_channel = setdiff ( 1:124, bad_channel );

    data_ch = data(good_channel, :);
    size_ch= size(data_ch, 1); 
    size_time = size(data_ch, 2);

    %redefine data
    EEG.data = data_ch;
    EEG.nbchan = size(data_ch, 1);
    EEG.nbchans = size(data_ch, 1);

    if subj_no == 1
       EGI_GSN_124(bad_channel) = [];
    end

    % edit channel location 
    EEG.chanlocs = EGI_GSN_124;

    % re-remove bad channels
    EEG = trimOutlier(EEG, 1e-4, 100, Inf, 0);
    EEG = eeg_checkset( EEG );

    bc_idx(subj_no).good = EEG.chanlocs;

    try 
      bc_idx(subj_no).bad = EEG.chaninfo.removedchans;
    catch
      warning('No channel removed.');
      bc_idx(subj_no).bad = [];
    end

    bc_idx(subj_no).name = fileName

end

% Save it as one of the two.

% save([outputDir 'bc_idx_C_original_new_25p'], 'bc_idx')
save([outputDir 'bc_idx_Ain_original_new_25p'], 'bc_idx')


%% Make a bad channel list & Remove bad channels 

clear

% change file path
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = [outputDir 'C_w_124_bp_original_new\']
outputDir2 = [outputDir 'C_w_124_bp_bc_original_new\']
outputDir3 = 'D:\HealthyBrainNetwork_above11\channel_location\';


load('D:\HealthyBrainNetwork_above11\channel_location\ch_128.mat')
load(['D:\HealthyBrainNetwork_above11\EEGDownload4\common_wb_124_bc_idx_25p_original.mat']);
load(['D:\HealthyBrainNetwork_above11\EEGDownload4\bc_idx_C_original_new_25p.mat'])


% load subj_list
subj_list = dir([outputDir1 '*_*']);

for subj_no = 1:length(subj_list)
    
    cd(outputDir1)
    
    % load the file
    input_fileName = char(subj_list(subj_no).name);
    fileName = erase(input_fileName, '.mat')
    load(input_fileName)

    bad_sym = [];
    bad_combined = [];
    bad_combined2 = [];
    bad_channel = [];
    bad_urchan = [];
    bad = [];
    
    data = EEG.data;
    Fs = EEG.srate;
    
    % remove bad channel 
    bad = bc_idx(subj_no).bad; 
    
    for i = 1:length(bad)
        bad_urchan(i) = bad(i).urchan;
    end
    
    bad_sym = ch_128(2, bad_urchan);
    bad_sym = bad_sym(bad_sym ~= 999);
    
    bad_combined = unique([bad_sym bad_urchan]);
    
    % Insert the bad common channel list
    bad_common = common_bc_sym;
    
    bad_combined2 = unique([bad_combined bad_common]);
    
    bad_channel = bad_combined2;
    ch_num = 124-numel(bad_channel);
    
    bc_idx(subj_no).bad_trim = bad_combined;
    bc_idx(subj_no).bad_all = bad_channel;
    bc_idx(subj_no).ch_num = ch_num;
    
    good_channel = setdiff ( 1:124, bad_channel );
    
    bc_idx(subj_no).good_trim = good_channel;
    
    data_ch = data(good_channel, :);
    size_ch= size(data_ch, 1); 
    size_time = size(data_ch, 2);
    
    %redefine data
    EEG.data = data_ch;
    EEG.nbchan = size(data_ch, 1);
    EEG.nbchans = size(data_ch, 1);

    chanlocs = EEG.chanlocs;
    chanlocs(bc_idx(subj_no).bad_all) = [];

    EEG.chanlocs = chanlocs;

    save_path_bc_re = [outputDir2 fileName '_bc_re']
    save(save_path_bc_re , '-v7.3', 'EEG', 'fileName')

    bc_idx_last(subj_no).good = EEG.chanlocs;

    try 
      bc_idx_last(subj_no).bad = bc_idx(subj_no).bad_all;
    catch
      warning('No channel removed.');
      bc_idx_last(subj_no).bad = [];
    end

    bc_idx_last(subj_no).name = fileName

    clear 'ch_labels*', 'bad*', 'bc_sym', 'good_channel'

end

% Save it as one of the two.
% save([outputDir 'bc_idx_C_original_new_25p_last_93'], 'bc_idx_last')
% save([outputDir 'bc_idx_Ain_original_new_25p_last_93'], 'bc_idx_last')



