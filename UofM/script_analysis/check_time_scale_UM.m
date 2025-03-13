% made by Youngjai Park, 2023.04.04.

% idx = 1;
% filter = 'alpha';

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

if ~exist([save_path '/check_time_resol/matfile_' filter], "dir")
    mkdir([save_path '/check_time_resol/matfile_' filter]);
end

% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
% load([save_path '/UM_info_AE.mat']);
% subjname = subjname_ae;
% labels = labels_ae;
load([save_path '/UM_info.mat']);
state = 'EC';
load([load_path '/preprocessing/slicing/band_' filter '/' subjname{idx} '_st_' state '_' labels{idx} '_band_' filter '_prepr.mat']);
time_resol_list = [2000, 1000, 500, 334, 250, 200:-2:2];
time_resol_list = flip(time_resol_list)/1000;

L = 100;                         % Overlap between windows (%)

for t = 1:size(time_resol_list,2)
    time_resol = time_resol_list(t); % Time resolution (s) = Time window size.
    H = double(data);

    % Input data must be time x channel
    % output_data is time window averaged data

    size_window = fs*time_resol; % unit: number of data points
    window_slide = size_window*L/100; % unit: number of data points
    time_rel_phase = size(H,1);

    mid_time_pt = (size_window-1)/2;

    time_pt = zeros( length(0:window_slide:time_rel_phase - size_window) , 1);
    output_data = zeros( length(0:window_slide:time_rel_phase - size_window) , size(H,2) );

    index = 1;

    if size_window ~= 1
        origin = load([save_path '/check_time_resol/matfile_' filter '/' subjname{idx} '_time_resolution_2ms.mat']);
        H = double(origin.output_data); % 나중에 어떤 데이터 쓸지 세팅하자.

        for i=0:window_slide:time_rel_phase - size_window

            time_pt(index) = mid_time_pt;

            output_data(index,:) = mean(H(i+1:size_window+i, :),1);

            mid_time_pt = mid_time_pt + window_slide;

            index = index + 1;

        end
    else
        output_data = cal_rel_phase_v3(H);
    	time_pt = 0:1:time_rel_phase-1;
    end
    save([save_path '/check_time_resol/matfile_' filter '/' subjname{idx} '_time_resolution_' num2str(time_resol*1000) 'ms.mat'], ...
        'output_data', 'time_pt', 'fs', 'window_slide', 'time_resol', 'chanlocs');
end