time_moving=1000;
time_window=2000;
fdat_data = chandata_real_60';
[Fdatm_phase, Fdatm_rel_p, Fdatm_rel_p_ave, Fdatm_phase_ave] = rel_phase_cal_mov( fdat_data, time_moving, time_window);

rel_p =  double(Fdatm_phase_ave);

k=22;

time_range_2 = [  24  26  40 44 46];
rel_p_mean_2 = mean(rel_p(time_range_2,:),1);

figure(2);
topoplot_general_test(rel_p_mean_2(:,:), chanloc_cord,'smooth',28);
title(  [ '\rm' num2str(time_now) 'S'; ] ,'fontsize',16 )
drawnow

time_range_1 = [ 140 231 233  239 241];
rel_p_mean_1 = mean(rel_p(time_range_1,:),1);

figure(1);
topoplot_general_test(rel_p_mean_1(:,:), chanloc_cord,'smooth',28);
title(  [ '\rm' num2str(time_now) 'S'; ] ,'fontsize',16 )
drawnow


% Frontoparietal:
% TR [6-17]: 8.4-33.6s (00:08 - 00:34)
% TR [26-28]: 50.4-56.7s (00:50 - 00:57)
% TR [78-83]: 159.6-172.2s (02:39 - 02:52)
% TR [97]: 199.5-201.6s (03:19 - 03:22)
% 9:2:33 51:2:55 161:2:171
% 8:34 50:57 159:172

time_range_4 = [ 8:2:34 50:2:57 ];
rel_p_mean_4 = mean(rel_p(time_range_4,:),1);

figure(4);
topoplot_general_test(rel_p_mean_4(:,:), chanloc_cord,'smooth',28);
%title(  [ '\rm' num2str(time_now) 'S'; ] ,'fontsize',16 )
drawnow

% Default:
% TR [8]: 12.6-14.7s (00:12 - 00:15)
% TR [12]:  21-23.1s (00:21 - 00:23)
% TR [19]: 35.7-37.8s (00:35 - 00:38)
% TR [25]: 48.3-50.4s (00:48 - 00:50)
% TR [38]: 75.6-77.7s (01:15 - 01:18)
% TR [41-42]: 81.9-86.1s (01:21 - 01:26)
% TR [51]: 102.9-105s (01:42 - 1:45)
% TR [61-62]: 123.9-128.1s (02:04 - 02:08)
% TR [67-68]: 136.5-140.7s (02:16 - 02:21)
% TR [78-80]: 159.6-165.9s (02:39 - 02:46) 
% TR [103]: 212.1-214.2s (03:32 - 03:34)
% TR [105]: 216.3-218.4s (03:36 - 03:38)

time_range_3 = [ 13:14 22 36:37 49 76:77 83:85 104 127 137:140 161:165 213 217];
rel_p_mean_3 = mean(rel_p(time_range_3,:),1);

figure(3);
topoplot_general_test(rel_p_mean_3(:,:), chanloc_cord,'smooth',25);
%title(  [ '\rm' num2str(time_now) 'S'; ] ,'fontsize',16 )
drawnow
 
