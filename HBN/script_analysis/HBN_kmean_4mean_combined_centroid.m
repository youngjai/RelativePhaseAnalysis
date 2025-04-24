
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load relatvie phase files 
% 2. Create all topo vector
% 3. Run K-means clusteirng 
% 4. Save centroids and clustered indices
%
% 2024. 02. 24. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all

cd 'D:\HealthyBrainNetwork_above11\EEGDownload4'

outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = [outputDir '\Ain_w_124_alpha_rp_eo_4mean_original_new'];

% make a folder
mkdir(outputDir1)


% change file path
file_path ='D:\HealthyBrainNetwork_above11\EEGDownload4\Ain_w_124_bp_bc_alpha_eo_original_new\';
file_list = dir([file_path '\*.mat']);

load([outputDir 'topo_ref_C_original.mat']) % minimu topo size reference

% must check the size of all topo_vector
topo_size2 = 31758;
time = 168; %ec 336 / eo 168 (*100 ms) the number of frame
all_topo_vector = zeros((time*length(file_list)), topo_size2);
topo_ref = topo_ref_C_original;

for subj_no = 1:length(file_list)
    cd(file_path) 
    input_fileName = char(file_list(subj_no).name);
    fileName = erase(input_fileName, '.mat') 

    load(input_fileName, 'H', 'time_moving', 'time_window', 'Fs', 'rel_p', 'rel_phase_w_mean','chan_coord_xy','time_all', 'time_moving', 'time_window', 'topo','smooth')
    
    shiftpreset = 0;
    % 0 if it is human
    % 1 if it is monkey
    
    % vectorize data
    topo_size = length(topo);
    
    temp1 = topo{1};
    topo_idx = isnan(temp1);
    [topo_idx_x,topo_idx_y] = find(topo_idx == 0);
    topo_vector = zeros(topo_size,topo_size2);
    
    % topo_vector: time point x frame vector
    for i=1:topo_size    
        T = topo{i};
        temp2  = T;
        temp2(topo_ref == 1) = NaN;
        T = temp2;
        T(isnan(T)) = [] ;
        topo_vector(i,:) = T;
    end
    
    all_topo_vector((time*(subj_no-1)+1):(time*subj_no), 1:size(topo_vector,2)) = topo_vector;

end


name = 'Ain_eo_original_new';
savepath = [outputDir 'all_topo\' name '_all_topo' ] 
save(savepath, 'all_topo_vector','name', '-v7.3');



fileName = name;

shiftpreset = 0;


%% K-mean clustering

load('D:\Dropbox\Moonbrainlab\HBN\centroid_topo_ref\comb_centroids_20240311/combined_centroids_20240311.mat');

K=4;
test_num = 10;

mask = zeros([K,length(topo_idx_comb)]);
for k = 1:K
    mask(k,:) = C_comb.EOEC(k,topo_idx_comb);
end

subj = 1;

tic;
all_topo_vector_new = all_topo_vector(:,topo_idx_HBN2comb);
toc;


tic;
[IDX,C,SUMD,D] = kmeans_with_mask_v2(mask, all_topo_vector_new);
toc;


% time windowing0
unit = 0.1;

% caculate mode, dwell time, transition
properties = cal_transition_prop_v2(IDX, K, unit)


save_path = [outputDir1 '\' fileName '_4mean_comb_cen_cal' ]
save( save_path, 'IDX', 'C', 'SUMD', 'D', 'properties');



