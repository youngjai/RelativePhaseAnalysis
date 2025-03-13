% made by youngjai 2023.03.07.
clear all; close all;

% sampling rate : 500 Hz
fs = 500;
% number of channels : 128
n_ch = 128;
 
load_path = [pwd() '/..']; % Project location
save_path = [load_path '/results_yjp'];

% 'event_info.mat' contains 'eventname', 'eventtime', 'subjname'
load([load_path '/matdata500/event_info.mat']);
logic = contains(subjname, 'UM');
subjname = subjname(logic);
eventtime = eventtime(logic,:);
n_sub = length(subjname);
% labels with matching subjname gender
% 1  UM_1  AE_64ch     M
% 2  UM_3  normal      M
% 3  UM_4  AE          M
% 4  UM_6  normal_64ch F
% 5  UM_7  AE          F
% 6  UM_8  BS          M
% 7  UM_9  AE          M
% 8  UM_10 normal      F
% 9  UM_11 normal      F
% 10 UM_12 BS          F
% 11 UM_13 BS          F
% 12 UM_14 BS          F
% 13 UM_16 normal      F
% 14 UM_18 BS          M * Caution!
% 15 UM_17 normal      F * Caution!
% 16 UM_20 normal      F
% 17 UM_21 BS          M
% 18 UM_22 normal      M (no EO state)!
% 19 UM_23 normal      M
% 20 UM_24 normal      M
labels = {'AE_64ch', 'normal', 'AE', 'normal_64ch', 'AE', 'BS', 'AE', ...
    'normal', 'normal', 'BS', 'BS', 'BS', 'normal', 'BS', ...
    'normal', 'normal', 'BS', 'normal', 'normal', 'normal'}';
gender = {'M', 'M', 'M', 'F', 'F', 'M', 'M', 'F', 'F', 'F', 'F', 'F', ...
    'F', 'M', 'F', 'F', 'M', 'M', 'M', 'M'}';

% Cautions! the time of events is not correct.
% You should check the "Note_McDonnell_AC_UC_TB_3_30_16.pdf" file.
% eventtime column:  CG1e(6), Dbg1(7), LOC(11), Dbg2(12), De2(13), ROC(15), CG2s(16)
%          3 UM_4 :  745000  '881500' '1211500' '1361500' '6434500' "7296000?" "7368500?"
%          5 UM_7 :  784665  1034077   1343217   1616151   6965335   7698928   7902945
%          6 UM_8 :  715167   836387   1064402   1279768   6507031   7735805   7854553
%          7 UM_9 :  780201  1062592   1334508  '1484508'  6807689   8341416   8506886
%         10 UM_12:  662728   943946   1199188  '1349188'  6750121   7375026   8201477
%         11 UM_13:  382511  1034941   1263237   1485780   6822565   7767088   7879535
%         12 UM_14:  749298  1217061   1574899   1669496   7005489   7984032   8063579
%         14 UM_18:  790275   994393   1235909   1458655   6840000  '7305334'  8001190
%         17 UM_21: 1014244  1167338   1496456  '1646456'  6894732   8739529   8806100
% ten values are different between the note file and the data file
% UM_inf_AE.mat use the time based on the note file except N/A
% "" time data is N/A in the pdf file so, we just use raw data

% UM_22 dosen't have EO time point. We just use 5 min data before ECs - 1sec.
eventtime(18,[1,2]) = [42215, 192215];

eventtime(3,[7,11,12,13]) = [881500, 1211500, 1361500, 6434500];
eventtime(7,12) = 1484508;
eventtime(10,12) = 1349188;
eventtime(14,15) = 7305334;
eventtime(17,12) = 1646456;

save([save_path '/UM_info.mat'], 'fs', 'subjname', 'eventtime', ...
    'eventname', 'n_sub', 'labels', 'n_ch', 'gender');

% save anesthesia subjects' information
n_sub_ae = 0;
for idx = 1:n_sub
    if ~contains(labels{idx}, {'normal', '64ch'})
        n_sub_ae = n_sub_ae+1;
    end
end

subjname_ae = cell(1,n_sub_ae);
eventtime_ae = zeros(n_sub_ae,size(eventtime,2));
labels_ae = cell(1,n_sub_ae);
gender_ae = cell(1,n_sub_ae);
idx_ae = 0;
for idx = 1:n_sub
    if ~contains(labels{idx}, {'normal', '64ch'})
        idx_ae = idx_ae+1;
        subjname_ae{idx_ae} = subjname{idx};
        eventtime_ae(idx_ae,:) = eventtime(idx,:);
        labels_ae{idx_ae} = labels{idx};
        gender_ae{idx_ae} = gender{idx};
    end
end

% start time and end time
slicetime_ae = zeros([n_sub_ae,6,2]);
cnt = 1;
% EO  : from eventtime_ae(idx, 1) to eventtime_ae(idx, 1)+5min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,1)+1;
slicetime_ae(:,cnt,2) = eventtime_ae(:,1)+5*60*fs;
cnt = cnt+1;
% EC  : from eventtime_ae(idx, 3) to eventtime_ae(idx, 3)+5min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,3)+1;
slicetime_ae(:,cnt,2) = eventtime_ae(:,3)+5*60*fs;
cnt = cnt+1;
% % AI1 : from eventtime_ae(idx, 7)+3min. to eventtime_ae(idx, 7)+8min.
% slicetime_ae(:,cnt,1) = eventtime_ae(:,7)+1+3*60*fs;
% slicetime_ae(:,cnt,2) = eventtime_ae(:,7)+8*60*fs;
% cnt = cnt+1;
% LOC : from eventtime_ae(idx,11) to eventtime_ae(idx,11)+3min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,11)+1;
slicetime_ae(:,cnt,2) = eventtime_ae(:,11)+3*60*fs;
cnt = cnt+1;
% % AI2 : from eventtime_ae(idx,12) to eventtime_ae(idx,12)+5min.
% slicetime_ae(:,cnt,1) = eventtime_ae(:,12)+1;
% slicetime_ae(:,cnt,2) = eventtime_ae(:,12)+5*60*fs;
% cnt = cnt+1;
% Anes : from eventtime_ae(idx,12)+10min to eventtime_ae(idx,12)+15min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,12)+1+10*60*fs;
slicetime_ae(:,cnt,2) = eventtime_ae(:,12)+15*60*fs;
cnt = cnt+1;
% BS  : refer to "timestemp_bs{2,idx}" and run "preprocessing_slicing_BS_UM.m"
% DS  : refer to "timestemp_bs{1,idx}" and run "preprocessing_slicing_BS_UM.m"
% DA  : from eventtime_ae(idx,13)-10min. to eventtime_ae(idx,13)-5min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,13)+1-10*60*fs;
slicetime_ae(:,cnt,2) = eventtime_ae(:,13)-5*60*fs;
cnt = cnt+1;
% ROC : from eventtime_ae(idx,20) to eventtime_ae(idx,20)+5min.
slicetime_ae(:,cnt,1) = eventtime_ae(:,20)+1;
slicetime_ae(:,cnt,2) = eventtime_ae(:,20)+5*60*fs;
    

filter_names = {'bf', 'delta', 'theta', 'alpha', 'low beta', 'high beta', ...
    'gamma', 'smr', 'alphab', 'thealp', 'lfbb', 'beta', ...
    '[1-2]', '[2-3]', '[3-4]', '[4-6]', '[6-8]', '[8-10]', '[10-12]', ...
    '[12-16]', '[16-20]', '[20-24]', '[24-28]', ...
    '[28-32]', '[32-40]', '[40-48]', '[1-12]'};
bands = [ 1 50; ... % [ 1-50] full band
          1  4; ... % [ 1- 4] delta band
          4  8; ... % [ 4- 8] theta band
          8 12; ... % [ 8-12] alpha band
         12 20; ... % [12-20] low beta band
         20 30; ... % [20-30] high beta band
         30 50; ... % [30-50] gamma band
         12 15; ... % [12-15] SMR (sensorimotor rhythm) band
          8 15; ... % [ 8-15] alpha band, broadly defined
          4 15; ... % [ 4-15] theta + alpha band
          1 15; ... % [ 1-15] low frequency broad band
         15 25; ... % [15-25] beta band
          1  2; ... % [ 1- 2] band
          2  3; ... % [ 2- 3] band
          3  4; ... % [ 3- 4] band
          4  6; ... % [ 4- 6] band
          6  8; ... % [ 6- 8] band
          8 10; ... % [ 8-10] band
         10 12; ... % [10-12] band
         12 16; ... % [12-16] band
         16 20; ... % [16-20] band
         20 24; ... % [20-24] band
         24 28; ... % [24-28] band
         28 32; ... % [28-32] band
         32 40; ... % [32-40] band
         40 48; ... % [40-48] band
          1 12];    % [ 1-12] band

save([save_path '/UM_info_AE.mat'], 'fs', 'subjname_ae', 'eventtime_ae', 'eventname', ...
    'n_sub_ae', 'labels_ae', 'n_ch', 'filter_names', 'bands', 'slicetime_ae', ...
    'gender_ae');

% save anesthesia subjects' information
n_sub_bs = 0;
for idx = 1:n_sub
    if matches(labels{idx}, 'BS')
        n_sub_bs = n_sub_bs+1;
    end
end

subjname_bs = cell(1,n_sub_bs);
eventtime_bs = zeros(n_sub_bs,size(eventtime,2));
labels_bs = cell(1,n_sub_bs);
gender_bs = cell(1,n_sub_bs);
idx_bs = 0;
for idx = 1:n_sub
    if matches(labels{idx}, 'BS')
        idx_bs = idx_bs+1;
        subjname_bs{idx_bs} = subjname{idx};
        eventtime_bs(idx_bs,:) = eventtime(idx,:);
        labels_bs{idx_bs} = labels{idx};
        gender_bs{idx_bs} = gender{idx};
    end
end

% Deep suppresion, and Burst suppression times
timestemp_bs = cell(2,n_sub_bs);
%         Deep suppresion 
% UM_8  : 1680000-1830000
% UM_12 : 2754000-2814000, 3084000-3174000
% UM_13 : 2400000-2550000
% UM_14 : 3060000-3150000, 3180000-3240000
% UM_18 : 2010000-2160000
% UM_21 : 2190000-2340000
timestemp_bs{1,1}(1:2) = [1680000, 1830000];
timestemp_bs{1,2}(1:2) = [2754000, 2814000];
timestemp_bs{1,2}(2,1:2) = [3084000, 3174000];
timestemp_bs{1,3}(1:2) = [2400000, 2550000];
timestemp_bs{1,4}(1:2) = [3060000, 3150000];
timestemp_bs{1,4}(2,1:2) = [3180000, 3240000];
timestemp_bs{1,5}(1:2) = [2010000, 2160000];
timestemp_bs{1,6}(1:2) = [2190000, 2340000];

% Burst suppresion 
% UM_8, UM_12, UM_13, UM_14, UM_18, UM_21
%    6,    10,    11,    12,    14,    17
%    3,     5,     6,     7,     8,     9
%    1,     2,     3,     4,     5,     6
% in the excel file, [start end] times are recorded
% UM_8, UM_12, UM_13, UM_14, UM_18, UM_21
for idx_bs = 1:n_sub_bs
    timestemp_bs{2,idx_bs} = xlsread([save_path '/ppt_materials/burst_time.xlsx'], subjname_bs{idx_bs});
end

save([save_path '/UM_info_BS.mat'], 'fs', 'subjname_bs', 'eventtime_bs', 'eventname', ...
    'n_sub_bs', 'labels_bs', 'timestemp_bs', 'n_ch', 'filter_names', 'bands', 'gender_bs');

% 128 channel subjects
idx_ch128 = ~contains(labels, '64ch');
subjname = subjname(idx_ch128);
n_sub = size(subjname,1);
eventtime = eventtime(idx_ch128,:);
gender = gender(idx_ch128);
labels = labels(idx_ch128);

save([save_path '/UM_info_ch128.mat'], 'fs', 'subjname', 'eventtime', ...
    'eventname', 'n_sub', 'labels', 'n_ch', 'gender', 'filter_names', 'bands', ...
    'slicetime_ae','timestemp_bs');

