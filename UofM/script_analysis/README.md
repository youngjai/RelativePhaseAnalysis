# README.md
% 2023.03.29.

## The order of run codes when we want to pre-process the data

1. making_info.m              % to make information about the data
% the script makes following variables
% fs : sampling rate
% subjname : 'UM_*'
% eventtime : when the task run
% eventname : the name series of the tasks
% n_sub : the number of subjcets
% labels : normal, AE(anesthesia), BS(burst suppression)

% subj  : UM_4, UM_7, UM_8, UM_9, UM_12, UM_13, UM_14, UM_18, UM_21
% idx   :    3,    5,    6,    7,    10,    11,    12,    14,    17
% idx_ae:    1,    2,    3,    4,     5,     6,     7,     8,     9
% idx_bs:    x,    x,    1,    x,     2,     3,     4,     5,     6

2. make_idx.m                     % to use when to run scripts
% three types of version (You should check the for loop statement!)
% * McDonnell UofM    subjects (n_sub   , eventtime   , labels      )
% * Deep anesthesia   subjects (n_sub_ae, eventtime_ae, labels_ae   )
% * Burst suppression subjects (n_sub_bs, eventtime_bs, timestemp_bs)

3. preprocessing_slicing_UM.m % to preprocess
  a) slice data with adding buffer of 2 sec
  b) band pass filter everything to detrend
  c) rereference (average referencing scheme) the data
  d) band pass filter the data for each power-band
  e) save the data
% EO  : from eventtime_ae(idx, 1) to eventtime_ae(idx, 1)+5min.
% EC  : from eventtime_ae(idx, 3) to eventtime_ae(idx, 3)+5min.
% AI1 : from eventtime_ae(idx, 7)+3min. to eventtime_ae(idx, 7)+8min.
% LOC : from eventtime_ae(idx,11) to eventtime_ae(idx,11)+3min.
% AI2 : from eventtime_ae(idx,12) to eventtime_ae(idx,12)+5min.
% BS  : refer to "timestemp_bs{2,idx}" and run "preprocessing_slicing_BS_UM.m"
% DS  : refer to "timestemp_bs{1,idx}" and run "preprocessing_slicing_BS_UM.m"
% DA  : from eventtime_ae(idx,13)-5min. to eventtime_ae(idx,13)   (UM_4,7,8,9,18,21)
% DA  : from eventtime_ae(idx,13)+5min. to eventtime_ae(idx,13)+10min. (UM_12,13,14)
% ROC : from eventtime_ae(idx,20) to eventtime_ae(idx,20)+5min.

4. integrating_data_UM.m      % to concatenate all slicing data
5. pspectrogram_UM.m          % to draw spectrogram to check primary


## The order of run codes when we want to measure dyanmics to check an optimal time window size
1. make_idx_filter.m          % to make a time series of relative phases according to a given time window (check_time_scale_UM.m)
                              % to calculate the correlation coefficient between the imaginary coherence and the relative phase (correlation_imc_relp_UM.m)
                              % to calculate the switching behavior of the correlation dynamics (calculate_imc_relp_dynamcis_UM.m)
                              % to calculate the ple matrix (calculate_ple_of_corr_UM.m)
                              % to calculate the interevent time distribution (draw_interevent_time_distribution.m)

2. make_idx_filter_measure.m  % to measure ApEn(Approximate Entropy), Shannon's entropy, Correlation Dimension, Standard deviation, Phase lag Entropy

3. make_idx_filter_1.m        % to calculate the mean of ple according to channel (calculate_avg_measure_of_rel_p_UM.m)


4. make_filter.m              % to draw a figure (draw_time_sacle_UM.m)
                              % to draw a figure (draw_time_scale_of_ple_UM.m)
5. make_idx.m                 % to draw a figure (draw_gif_time_sacle_UM.m)
                              % to draw a figure (draw_gif_time_sacle_UM_1.m)
