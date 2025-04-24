%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load EEG files
% 2. Bandpass filtering - alpha
% 3. Seperate by event time
% 4. Calculate relative phase

% 2024. 01. 16. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% Bandpass filtering
% change file path
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = [outputDir 'C_w_124_bp_bc_original_new\']
outputDir2 = [outputDir 'C_w_124_bp_bc_alpha_original_new\']
outputDir3 = [outputDir 'C_w_124_bp_bc_alpha_event_original_new\']
mkdir([outputDir2])
mkdir([outputDir3])

% load subj_list
subj_list = dir([outputDir1 '*_*']);

% Banspass filtering for the specific frequency ranges
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

    % bandpass fitering (alpha)
    EEG = pop_eegfiltnew(EEG, 'locutoff',8,'hicutoff',12,'plotfreqz',1);
    close gcf

    % event extract
    data = EEG.data;
    Fs = EEG.srate;

   % event time 
    event = EEG.event;
    event_sample = [event.sample];
    event_type = transpose({event(:).type});

%     event(1) = [];
    
    event_time = event_sample*Fs/60;
    
    RS_total.event(subj_no).type = event_type;
    RS_total.event(subj_no).sample = event_sample;
    RS_total.event(subj_no).time = event_time;
    
    % data epoch
    RS_total.data(subj_no).epoch = {};
    RS_total.data(subj_no).epoch_fix = {};
    RS_total.data(subj_no).EO = {};
    RS_total.data(subj_no).EC = {}; 
    RS_total.state(subj_no).EO= []; 
    RS_total.state(subj_no).EC= [];  
    RS_total.state(subj_no).EO_event = []; 
    RS_total.state(subj_no).EC_event = []; 
    
    event_sample_temp = event_sample;
    
    for es = 1:11
        event_sample_new(es) = event_sample_temp(es)-(event_sample_temp(1)-1);
    end
    
    event_sample_new = event_sample_new';
    event_sample = event_sample_new;
    
    tic
    for ch = 1:size(data, 1)
        for es = 1:11 %event_sample
    
            if es < 11
                temp = data(ch,(event_sample(es):event_sample((es)+1))); 
                temp = nonzeros(temp);
                RS_total.data(subj_no).epoch(ch).dat(es, 1:length(temp)) = temp;
            elseif es == 11
                continue
            end %if
    
        end % es
    
        for es = 1:10 
            temp = RS_total.data(subj_no).epoch(ch).dat(es, 1:20001);
            RS_total.data(subj_no).epoch_fix(ch).dat(es, 1:length(temp)) = temp;
        end 
    
         % data epoch for state
        for en = 1:10 
            
            if ch == 1 % event time record
            
                if floor(en/2) == en/2 % even number == EC
                   RS_total.data(subj_no).EC(ch).dat(en/2, :) = RS_total.data(subj_no).epoch_fix(ch).dat((en), 1:20001);
                   temp2 = transpose(nonzeros(RS_total.data(subj_no).EC(ch).dat(en/2, :)));
    
    
                   if en/2 == 1
                        RS_total.state(subj_no).EC(ch, 1:length(temp2)) = temp2;
                        RS_total.state(subj_no).EC_event(en/2, :) = length(temp2);   
                   else
                        temp3 = length(nonzeros(RS_total.state(subj_no).EC(ch, :)));
                        RS_total.state(subj_no).EC(ch, (temp3+1):temp3+length(temp2)) = temp2;
                        RS_total.state(subj_no).EC_event(en/2, :) = length(RS_total.state(subj_no).EC(ch, :)); 
                   end
    
                else % odd number == EO
                   RS_total.data(subj_no).EO(ch).dat((en+1)/2, :) = RS_total.data(subj_no).epoch_fix(ch).dat((en), 1:10001);
                   temp2 = transpose(nonzeros(RS_total.data(subj_no).EO(ch).dat((en+1)/2, :)));
    
                   if (en+1)/2 == 1
                        RS_total.state(subj_no).EO(ch, 1:length(temp2)) = temp2;
                        RS_total.state(subj_no).EO_event((en+1)/2, :) = length(temp2);  
                   else
                        temp3 = length(nonzeros(RS_total.state(subj_no).EO(ch, :)));
                        RS_total.state(subj_no).EO(ch, (temp3+1):temp3+length(temp2)) = temp2;
                        RS_total.state(subj_no).EO_event((en+1)/2, :) = length(RS_total.state(subj_no).EO(ch, :));
                   end
    
                end  
    
        
            else 
                if floor(en/2) == en/2 % even number == EC
                   RS_total.data(subj_no).EC(ch).dat(en/2, :) = RS_total.data(subj_no).epoch_fix(ch).dat((en), 1:20001);
                   temp2 = transpose(nonzeros(RS_total.data(subj_no).EC(ch).dat(en/2, :)));
    
    
                   if en/2 == 1
                        RS_total.state(subj_no).EC(ch, 1:length(temp2)) = temp2; 
                   else
                        temp3 = length(nonzeros(RS_total.state(subj_no).EC(ch, :)));
                        RS_total.state(subj_no).EC(ch, (temp3+1):temp3+length(temp2)) = temp2;
                   end
    
                else % odd number == EO
                   RS_total.data(subj_no).EO(ch).dat((en+1)/2, :) = RS_total.data(subj_no).epoch_fix(ch).dat((en), 1:10001);
                   temp2 = transpose(nonzeros(RS_total.data(subj_no).EO(ch).dat((en+1)/2, :)));
    
                   if (en+1)/2 == 1
                        RS_total.state(subj_no).EO(ch, 1:length(temp2)) = temp2;
                   else
                        temp3 = length(nonzeros(RS_total.state(subj_no).EO(ch, :)));
                        RS_total.state(subj_no).EO(ch, (temp3+1):temp3+length(temp2)) = temp2;
                   end
    
                end  
    
            end % event time record
        end % en
    end % ch 

    save_path_alpha = [outputDir2 fileName '_alpha']
    save(save_path_alpha, '-v7.3', 'EEG', 'fileName')

    ch_num = length(RS_total.data(subj_no).EO);


    % Change the segment you want to save.
    for i = 1:5
        for j = 1:ch_num
              EC(i).dat(j, :) = RS_total.data(subj_no).EC(j).dat(i, 2701:19500);
              EO(i).dat(j, :) = RS_total.data(subj_no).EO(j).dat(i, 1101:9500);
        end
    end

    save_path_alpha_event = [outputDir3 fileName '_event']
    save(save_path_alpha_event , '-v7.3', 'EO', 'EC', 'fileName')

    clear('EC', 'EO')
end


%% Relative phase 
clear

% change file path
outputDir = 'D:\HealthyBrainNetwork_above11\EEGDownload4\';
outputDir1 = [outputDir 'C_w_124_bp_bc_alpha_event_original_new\']
outputDir2 = [outputDir 'C_w_124_bp_bc_alpha_eo_original_new\']
outputDir3 = [outputDir 'C_w_124_bp_bc_alpha_ec_original_new\']
mkdir([outputDir2])
mkdir([outputDir3])

%load('D:\HealthyBrainNetwork_above11\EEGDownload4\bc_idx_C_original_new_25p_last_93.mat');
load('D:\HealthyBrainNetwork_above11\EEGDownload4\bc_idx_C_original_new_25p_last.mat');

bc_idx = bc_idx_last;


% load subj_list
file_list = dir([outputDir1 '*_*']);

smooth = 20; % spatial smoothness
tend_temp = 1; % portion of the time which you want to make the movie with
tend = 1./tend_temp;

EO_NAME = {'EO1', 'EO2', 'EO3', 'EO4', 'EO5'};
EC_NAME = {'EC1', 'EC2', 'EC3', 'EC4', 'EC5'};

for subj_no = 1:length(file_list)
   
    cd (outputDir1)
    
    input_fileName = char(file_list(subj_no).name);
    fileName = erase(input_fileName, '.mat') 

    load(input_fileName)

    load('D:\HealthyBrainNetwork_above11\channel_location\chan_coord_info.mat')

    tic
    for eo = 1:5
    
        H = EO(eo).dat(:, :)';
        Fs = 500;

        
        for i = 1:length(bc_idx(subj_no).good)
            ch_idx(i) = bc_idx(subj_no).good(i).urchan;
        end 
%         ch_idx = ch_idx - 3;

        chan_coord = chan_coord_org(ch_idx, :);
        
        chan_coord_xy = chan_coord;

        % Compute relative phase
        [rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta] = cal_rel_phase_v2(H) ;
        
        % Moving time window
        time_moving = 50;  %  number of time samples, time_moving / Fs = moving time in seconds
        time_window = 50;  %  number of time samples, time_window / Fs = window time in seconds
        dt=1;
        %Fs = 500;
        
        [ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);
        
        rel_p = double( rel_phase_w_mean(1:round(size(rel_phase_w_mean,1)./tend),: ) ) ;
        
        
        %vidfile = VideoWriter([ file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_relp' ] ,'MPEG-4');
        %vidfile.FrameRate = round(Fs/time_moving/dt);
        %vidfile.Quality = 75;
        %open(vidfile);    
        %FrameRate = round(Fs/time_moving/dt);
        
        time_dt = 1/Fs*time_moving;
        time_now = 0 + (time_window)/2 /Fs;
        time_all = [1:dt:round(size(rel_p,1)./tend )]';
        
        
        formatSpec = '%.2f';
        topo = cell(size(time_all,1),1);
        
        tic
        for i=1:dt:(size(rel_p,1))
            
        figure(1);
        topo{i} = topoplot_general_test(rel_p(i,:)', chan_coord_xy(:,1:2),'smooth',smooth,'scatter',1);
        title(  [ '\rm' num2str(time_now,formatSpec) ' S'; ] ,'fontsize',16 )
        drawnow
        time_all(i)=time_now;
        F(i) = getframe(gcf); 
        %writeVideo(vidfile,F(i));
        time_now = time_now + dt*time_dt;
        end
        %close(vidfile)
        toc
        
        save_path_rp = [outputDir2 fileName '_' EO_NAME{eo} '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth)  '_rp']
        save(save_path_rp, 'H', 'Fs', 'rel_p', 'rel_phase_w_mean', 'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','-v7.3' );
        
        ch_idx = [];
            
    end % eo(1:5)

    toc
    tic
    for ec = 1:5
    
        H = EC(ec).dat(:, :)';
        Fs = 500;


        for i = 1:length(bc_idx(subj_no).good)
            ch_idx(i) = bc_idx(subj_no).good(i).urchan;
        end 
%         ch_idx = ch_idx - 3;
        
        chan_coord = chan_coord_org(ch_idx, :);
        
        chan_coord_xy = chan_coord;

        % Compute relative phase
        [rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta] = cal_rel_phase_v2(H) ;
        
        % Moving time window
        time_moving = 50;  %  number of time samples, time_moving / Fs = moving time in seconds
        time_window = 50;  %  number of time samples, time_window / Fs = window time in seconds
        dt=1;
        %Fs = 500;
        
        [ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);
        
        rel_p = double( rel_phase_w_mean(1:round(size(rel_phase_w_mean,1)./tend),: ) ) ;
        
        
        %vidfile = VideoWriter([ file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_relp' ] ,'MPEG-4');
        %vidfile.FrameRate = round(Fs/time_moving/dt);
        %vidfile.Quality = 75;
        %open(vidfile);    
        %FrameRate = round(Fs/time_moving/dt);
        
        time_dt = 1/Fs*time_moving;
        time_now = 0 + (time_window)/2 /Fs;
        time_all = [1:dt:round(size(rel_p,1)./tend )]';
        
        
        formatSpec = '%.2f';
        topo = cell(size(time_all,1),1);
        
        tic
        for i=1:dt:(size(rel_p,1))
            
        figure(1);
        topo{i} = topoplot_general_test(rel_p(i,:)', chan_coord_xy(:,1:2),'smooth',smooth,'scatter',1);
        title(  [ '\rm' num2str(time_now,formatSpec) ' S'; ] ,'fontsize',16 )
        drawnow
        time_all(i)=time_now;
        F(i) = getframe(gcf); 
        %writeVideo(vidfile,F(i));
        time_now = time_now + dt*time_dt;
        end
        %close(vidfile)
        toc
        
        
        save_path_rp = [outputDir3 fileName '_' EC_NAME{ec} '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth)  '_rp']
        save(save_path_rp, 'H', 'Fs', 'rel_p', 'rel_phase_w_mean', 'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','-v7.3' );
         
        ch_idx = [];
            
    end % ec(1:5)
    toc

end % subject

