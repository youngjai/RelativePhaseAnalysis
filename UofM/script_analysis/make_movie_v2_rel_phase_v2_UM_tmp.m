% clear all;

smooth = 20; % spatial smoothness
tend_temp = 1; % portion of the time which you want to make the movie with
tend = 1./tend_temp;

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

% 'event_info.mat' contains 'eventname', 'eventtime', and 'subjname'
% load([save_path '/UM_info_AE.mat']);
% subjname = subjname_ae;
% labels = labels_ae;
load([save_path '/UM_info_ch128.mat']);

% idx = 4;
% filter = 'alpha';
% filter = 4;
band_name = ['band_[' num2str(bands(filter,1)) '-' ...
    num2str(bands(filter,2)) ']'];
% state = 'EO';
state_list = {'EO', 'EC', 'LOC', 'AN', 'BS', 'DS', 'DA', 'ROC'};
st = find(matches(state_list, state));

% Define file name and load the file.
file_name = [subjname{idx} '_st_' state '_' band_name];
load([load_path '/preprocessing/completed_240715/' band_name '/' subjname{idx} ...
    '_wo_badchan_and_wo_reref_and_filtered.mat']);
data = data_state{st};

H = double(data);

chanlocs = chanlocs{st};
chan_coord_xy = [ -[chanlocs.Y]' [chanlocs.X]']; % xy should be changed due to the different definition

% Compute relative phase
[rel_phase_w, ~, ~, ~, ~, ~, ~] = cal_rel_phase_v3(H) ;

% Moving time window
time_moving = 10;  %  number of time samples, time_moving / fs = moving time in seconds
time_window = 10;  %  number of time samples, time_window / fs = window time in seconds
dt=1;

[ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);

rel_p = double( rel_phase_w_mean(1:round(size(rel_phase_w_mean,1)./tend),: ) ) ;

if ~exist([save_path '/movie_rel_phase'], "dir")
    mkdir([save_path '/movie_rel_phase']);
end
save_file_path = [save_path '/movie_rel_phase/' band_name '_20240723'];
% save_file_path = [save_path '/movie_rel_phase/' band_name '_20240402_wo_reref/20ms'];
if ~exist(save_file_path, "dir")
    mkdir(save_file_path);
end
save_file_path = [save_path '/movie_rel_phase/' band_name '_20240723/' num2str(time_moving/fs*1000)  'ms'];
if ~exist(save_file_path, "dir")
    mkdir(save_file_path);
end
if ~exist([save_file_path '/' subjname{idx}], "dir")
    mkdir([save_file_path '/' subjname{idx}]);
end

% %vidfile = VideoWriter([ file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_relp' ] ,'MPEG-4');
% vidfile = VideoWriter([save_path '/movie_rel_phase/band_' filter '_whole/' subjname{idx} '/' ...
%     file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_relp' ]);
% vidfile.FrameRate = round(fs/time_moving/dt);
% vidfile.Quality = 75;
% open(vidfile);    
% FrameRate = round(fs/time_moving/dt);

time_dt = 1/fs*time_moving;
time_now = 0 + (time_window)/2 /fs;
time_all = [1:dt:round(size(rel_p,1)./tend )]';


formatSpec = '%.3f';
topo = cell(size(time_all,1),1);

tic;
for i=1:dt:size(rel_p,1)
    disp(i);
    % figure(1);
%     topo{i} = topoplot_general_test(rel_p(i,:)', chan_coord_xy(:,1:2),'smooth',smooth,'scatter',1);
    % when you want to make movie, use topoplot_figure(topo{i}, chan_coord_xy(:,1:2), 'scatter',1);
    topo{i} = topoplot_general_test_without_figure(rel_p(i,:)', chan_coord_xy(:,1:2),'smooth',smooth,'scatter',1);
%     title(  [ '\rm' num2str(time_now,formatSpec) ' s'; ] ,'fontsize',16 );
%     set(gcf, 'Visible', 'off');
    % drawnow
    time_all(i)=time_now;
%     F(i) = getframe(gcf); 
%     writeVideo(vidfile,F(i));
    time_now = time_now + dt*time_dt;
end
% close(vidfile)
toc;

save( [save_file_path '/' subjname{idx} '/' ...
    file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth)  '_topo' ], ...
    'H', 'fs', 'rel_p', 'rel_phase_w_mean', 'chan_coord_xy','time_all', 'time_pt', ...
    'time_moving', 'time_window', 'topo','smooth','-v7.3' );

% save topo_vector
save_file_path = [save_file_path '/topo_vector'];
if ~exist(save_file_path, "dir")
    mkdir(save_file_path);
end

load('references/comb_centroids_20240311/combined_centroids_20240311.mat', 'topo_idx_comb');
topo_vector = zeros([size(rel_p,1), length(topo_idx_comb)]);
for t = 1:size(rel_p,1)
    topo_vector(t,:) = topo{t}(topo_idx_comb);
end
save([save_file_path '/topo_vector_st_' state '_' subjname{idx} '.mat'], ...
    'topo_vector', '-v7.3');