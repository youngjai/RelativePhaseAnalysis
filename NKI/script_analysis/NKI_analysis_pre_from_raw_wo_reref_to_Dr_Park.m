
%%
% NKI EEG-fMRI data analysis
%
%
% Hyoungkyu Kim 2023.05.18

% number of states in each subject
% 12 12 12 12 23 23 23 / 24 12 24 24 24 / 24 12 23 23 30 / 12 24 24 24 24
% 1:7 8:12 13:17 18:22


%%
% GA PA removal

clear all
close all
clc
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% loading file names in data directory
datapath=['E:\NKI_download\raw_data\'];
datafold=dir([datapath 'sub-*']);
tic
real_time=clock;
[num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm']
bandname={'d', 't', 'a1', 'a2', 'b', 'r', 'w', '10'};
bandname2={'DELTA', 'THETA', 'ALPHA', 'ALPHA_narrow', 'BETA', 'GAMMA', 'WHOLE', '10Hz', 'ALPHA', 'Wave10'};
band_freq = [0.5 4 ; 4 8 ; 8 12 ; 9 11 ; 15 30 ; 30 50; 1 50; 9.5 10.5; 8 13];   % spectrogram paper 0.05~50Hz
electrodeExclude = 32; % ECG channel
%%%%%%%%%%%%%%%%%%%%
band_no=3;  
%%%%%%%%%%%%%%%%%%%%
sf=5000;  % sampling frequency = 500 Hz
% width=5;  % for wavelet filter parameter
% Nsurro=20;
% tr_drop=[0 1 2];  % 0 : no trial drop, 1 : light trial drop, 2 : heavy channel drop 
ch_no_F=sort([1 2 60 38 39 46 3 32 17 33 4 47]); % Fp1 Fp2 Fp AF3 AF4 F5 F3 F1 Fz F2 F4 F6 
ch_no_P=sort([50 7 36 19 37 8 51]);    % P5 P3 P1 Pz P2 P4 P6 
ch_no_O=sort([44 31 45]);              % PO3 POz Po4  //  9 10 20 - O1 O2 Oz 
ch_Cz=sort([34 18 35 23 61 24]);       % C1 Cz C2 CP1 CPz CP2

reference_name={'Single', 'Average', 'Laplacian', 'Mastoid'};  %, 'AverageTB', 'LaplacianTB'};
ref_type_no=2;  % 1,2,3,4
ANCFlag=1; % 1 or 0
outputfolder_eeg = ['E:\NKI_download\preproc_self\eeg\'];
mkdir(outputfolder_eeg)
mkdir([outputfolder_eeg '\figures'])

No_Pz = [];   pz_idx = [];   poz_idx = [];
for sub_no=1 % 8:length(datafold) %  1:7 8:12 13:17 18:22
    ses_list = dir([datapath datafold(sub_no).name '\ses-*']);

    for ses_no=2 % 1:length(ses_list)

        st_list = dir([datapath datafold(sub_no).name '\' ses_list(ses_no).name '\eeg\*.set']);

        for st_no=9:length(st_list)
            filename = st_list(st_no).name;
            filepath = [datapath datafold(sub_no).name '\' ses_list(ses_no).name '\eeg\'];
            
            %%%%%%%%%%%%%%%%%%%%% name of the current file to use
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            EEG = pop_loadset('filename',filename,'filepath',filepath);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

            output_fileName = filename;
            output_fileName = strrep(output_fileName,'.set','');

            tmp_name=st_list(st_no).name;
            if contains(tmp_name,'off') == 1 % scanner OFF
                processing_idx = 0;
            elseif contains(tmp_name,'out') == 1
                processing_idx = -1; % outside of scanner
            else 
                processing_idx = 1;
            end
    
            %% Gradient Artifact Removal
%             EEG.data = double(EEG.data);
            
            % run Gradient artifact removal only in ON condition
            if processing_idx==1 %|| (processing_idx==0)
        
                if (ANCFlag == 0)
                    output_fileSET = fullfile(outputfolder_eeg,[output_fileName,'_01_gradient_noANC','.set']);
                    EEG = pop_fmrib_fastr(EEG,[],[],[],'R128',1,0,[],[],[],[],electrodeExclude,'0');
                    [~,EEG,~] = pop_newset([], EEG, 1,'setname',[EEG.setname,' | GA Removed (no ANC)'],'savenew',output_fileSET,'gui','off');

                elseif (ANCFlag == 1)
                    output_fileSET = fullfile(outputfolder_eeg,[output_fileName,'_gradient_ANC','.set']);
                    EEG = pop_fmrib_fastr(EEG,[],[],[],'R128',1,1,[],[],[],[],electrodeExclude,'0');  % auto -> 0
%                     [~,EEG,~] = pop_newset([], EEG, 1,'setname',[EEG.setname,' | GA Removed (ANC)'],'savenew',output_fileSET,'gui','off');
                end
            end % if processing_idx ==1
            
            ecgchan = 32;

            % divide files into two folders (but the same name '_PArm')
            if ANCFlag==0
                output_fileSET = fullfile(outputDir0,[input_fileName(1:length(input_fileName)-5),'_PArm','.set']);
            elseif ANCFlag==1
                output_fileSET = fullfile(outputfolder_eeg,[output_fileName,'_GAPArm','.set']);
            end

            EEG = pop_fmrib_qrsdetect(EEG,ecgchan,'qrs','no');
            EEG = pop_fmrib_pas(EEG,'qrs','median');
            [~,EEG,~] = pop_newset([], EEG, 1,'setname',[EEG.setname,' | GA PA Removed'],'savenew',output_fileSET,'gui','off');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %             %%%%%%% remove ECG, EOGL, EOGU - index 32 63 64
% %                 EEG = pop_rejchan(EEG, 'elec',[32 63 64],'threshold',0.001,'norm','off','measure','kurt');
% %                 [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% 
%             %%%%%%re-referencing
%                 no_of_channels=EEG.nbchan;
%                 if ref_type_no==2
%                     EEG = pop_reref( EEG, 1: no_of_channels, 'keepref', 'on' );
%                     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
%                 elseif ref_type_no==3
%                     for ch_location=1:no_of_channels
%                         XX(ch_location,1)=EEG.chanlocs(ch_location).X;
%                         YY(ch_location,1)=EEG.chanlocs(ch_location).Y;
%                         ZZ(ch_location,1)=EEG.chanlocs(ch_location).Z;
%                     end
%                     [Laplacian_ref,G,H] = laplacian_perrinX(EEG.data,XX,YY,ZZ);  
%                     EEG.data=Laplacian_ref;
%                     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
%                 elseif ref_type_no==4
%                     tmpx=EEG.data;
%                     for i=1:size(EEG.data,1)
%                         tmpx(i,:)=tmpx(i,:)-mastoid;
%                     end
%                     EEG.data=tmpx;
%                     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
%                 end
% 
%             % Spectrogram
%             pz_ox=0;
%             for ch_no=1:length(EEG.chanlocs)
%                 ch_name = EEG.chanlocs(ch_no).labels;
%                 if ch_name(1) == 'P' && ch_name(2) == 'z'
%                     pz_idx(sub_no, ses_no, st_no) = ch_no;
%                     pz_ox=1;
%                 end % if 
%             end % ch_no
% 
%             poz_ox=0;
%             if pz_ox == 0
%                 No_Pz(sub_no, ses_no, st_no) = 1;
%                 for ch_no=1:length(EEG.chanlocs)
%                     ch_name = EEG.chanlocs(ch_no).labels;
%                     if ch_name(1) == 'P' && ch_name(2) == 'O'
%                         poz_idx(sub_no, ses_no, st_no) = ch_no;
%                         poz_ox=1;
%                     end % if 
%                 end % ch_no
%             end % if
% 
%             plot_ch_idx = 0;
%             if pz_ox == 1
%                 plot_ch_idx = pz_idx(sub_no, ses_no, st_no); % [20 9 10 31 59 60];
%             elseif poz_ox == 1
%                 plot_ch_idx = poz_idx(sub_no, ses_no, st_no);
%             end 
% 
%             if plot_ch_idx>0
%                 tmp_data=double(EEG.data);
%                 
%                 time_min=1; % 100*sf
%                 time_max=length(tmp_data); % 192*sf;
%                 tmp_data(isnan(tmp_data))=0;
%                 NFFT=100000/20; sf=5000/20; window=10000/20; overlap=round(4*window/5); 
%                 
%                 [b, F ,T, P]=spectrogram( tmp_data(plot_ch_idx,time_min:time_max),window,overlap,NFFT,sf);
%                 
%                 figure1_a=figure;
%                 set (figure1_a, 'resize','off')
%                 % set (figure1_a,'DefaultTextFontName','Arial','DefaultTextFontSize',10, 'DefaultTextFontWeight','Normal');
%                 % set (figure1_a, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontWeight','Normal');
%                 % set (gca, 'FontName', 'Arial', 'FontSize', 12);
%                 set (figure1_a, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                             %         pos = [0 0 20 10];
%                 pos = [20 10 25 15];  % [0 0 40 24];
%                 set(figure1_a, 'Position', pos, 'Units', 'centimeters');
%                 % set(figure1_a, 'PaperPositionMode', 'auto')     
%         
%                 imagesc( smooth2a(log(abs(double(b(1:901,:)))),1), [-10 10]), axis xy, colormap(jet)
%                 
%                     set(gca,'FontSize',12,'YTick',[10 101 201 241 301 401 501 601 901],'YTickLabel', {'0.5','5','10','12', '15','20','25','30','45'});
%                     set(gca,'FontSize',12,'XTick',[1:50:size(b,2)],'XTickLabel', [0:20:400]); % {'1','5','10','15','20','25','30','35'}
%                     ylabel('Frequency', 'fontsize', 14)
%                     xlabel('Time (sec)', 'fontsize', 14)
%                     new_name = strrep([filename(1:end-4)],'_',' ');
%                     title([new_name], 'fontsize', 14)
%         
%                 print(figure1_a,'-dtiffn', '-r300', [outputfolder '\figures\' new_name])
%                 close gcf 
%                 clear T F
%             end % if
% 
%             %%%%%%%%%%%%%%%%%%%%% FILTERING STEP
% %                 band_no=4;
%                     EEG = pop_firws(EEG, 'fcutoff', band_freq(band_no,1), 'ftype', 'highpass', 'wtype', 'blackman', 'forder', 2750, 'minphase', 0);
%                     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
%                     EEG = pop_firws(EEG, 'fcutoff', band_freq(band_no,2), 'ftype', 'lowpass', 'wtype', 'blackman', 'forder', 2750, 'minphase', 0);
%                     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
%                     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% 
%             data=double(EEG.data)';
% 
%             tr_mk=[];
%             mk_count=0;
%             for mk_no=1:length(EEG.event)
%                 marker_tmp = EEG.event(mk_no).type;
%                 if  marker_tmp(1) == 'Q'
%                     mk_count = mk_count +1;
%                     tr_mk(mk_count,1) = EEG.event(mk_no).latency;
%                 end % if
%             end % mk_no
% 
%             smooth = 20;
%             tend = 1; % 1/3; % portion of the time which you want to make the movie with
%                 
%             % Define file name and load the file.
%             H = data(round(tr_mk(1)+1):round(tr_mk(end-1)),:);
%             for ch_no=1:length(EEG.chanlocs)
%                 chan_coord(ch_no,1)=EEG.chanlocs(ch_no).X;
%                 chan_coord(ch_no,2)=EEG.chanlocs(ch_no).Y;
%                 chan_coord(ch_no,3)=EEG.chanlocs(ch_no).Z;
%             end % ch_no
%             chan_coord_yx = [ chan_coord(:,2) chan_coord(:,1) chan_coord(:,3) ];
%                 
%             % Compute relative phase
%             [rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta] = cal_rel_phase_v2(H) ;
%                 
%             % Moving time window
%             time_moving = 25; % 100ms for 250Hz SF
%             time_window = 25; % 1050; % 50;
%             dt=1;
%             Fs = 250; % 5000;
%                 
%             [ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);
%                 
%             rel_p =  double(rel_phase_w_mean);
% 
%             % making topo and movie
%             vidfile = VideoWriter([ outputfolder '\figures\' new_name '_tm' num2str(time_moving) '_tw' num2str(time_window) ... 
%                                   '_sm' num2str(smooth) '_' num2str(band_freq(band_no,1)) 'to' num2str(band_freq(band_no,2)) 'Hz'] ,'MPEG-4');
%             vidfile.FrameRate = round(Fs/time_moving/dt); % 0.48 for 2100ms
%             vidfile.Quality = 75;
%             open(vidfile);    
%             FrameRate = round(Fs/time_moving/dt);
%                 
%             time_dt = 1/Fs*time_moving;
%             time_now = 0 + (time_window)/2/Fs;
%             time_all = [1:dt:round(size(rel_p,1)./tend )]';
%                 
%                 
%             formatSpec = '%.2f';
%             time_set = size(time_all,1); % 100;
%             topo = cell(time_set,1);  % cell(size(time_all,1),1);
%                 
%             tic
%             warning('off','all')
%             for i=1:dt:round(size(rel_p,1)/tend)
%                     
%                 fig_H=figure;
%                 topo{i} = topoplot_general_test(rel_p(i,:)', chan_coord_yx(:,1:2),'smooth',smooth,'scatter',1);
%                 title(  [ '\rm' num2str(time_now,formatSpec) ' S'; ] ,'fontsize',16 )
%                 drawnow
%                 time_all(i)=time_now;
%                 F(i) = getframe(gcf); 
%                 writeVideo(vidfile,F(i));
%                 time_now = time_now + dt*time_dt;
%                 close
%             end
%             close(vidfile)
% %             warning('on','all')
%             toc
% 
%             shiftpreset = 0;
%             % 1 if it is monkey
%             % 0 if it is human
%                 
%             % vectorize data
%             topo_size = length(topo);
%                 
%             temp1 = topo{1};
%             topo_idx = isnan(temp1);
%             [topo_idx_x,topo_idx_y] = find(topo_idx == 0);
%             topo_vector = zeros(topo_size,length(topo_idx_x));
%                 
%             % topo_vector: time point x frame vector
%             for i=1:topo_size
%                 T = topo{i};
%                 T(isnan(T)) = [] ;
%                 topo_vector(i,:) = T;
%             end
%                 
%                 
%             %% This step is to test for the optimal K 
%             % This step takes very long time, skip if you are certain about the number K that you want to try
%             tic
%             eva1 = evalclusters(topo_vector(end-200:end-1,:),'kmeans','gap','KList',[1:20]);
% %             eva2 = evalclusters(topo_vector(end-600:end-401,:),'kmeans','gap','KList',[1:10]);
% %             eva = evalclusters(rand(10000,1),'kmeans','gap','KList',[1:10]);
%             toc
%                 
%             %% K-mean clustering
%             K=eva1.OptimalK; % round((eva1.OptimalK+eva2.OptimalK)/2);
%                 
%             % k-mean clustering core
%                 
%                 
%             tic
%             disp( [ 'k=' num2str(i) ] )
%             [IDX,C,SUMD,D]=kmeans(topo_vector,K, 'distance', 'sqeuclidean','Replicates',5,'Display','final','emptyaction','drop');   
%                 
%             % [s,h] = silhouette(topo_vector,IDX);
%             [sil] = silhouette(topo_vector,IDX);
%             SIL = mean(sil);
%             toc
%                 
% % %           set(0, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition',  [-0.5 0.25 0.35 0.5])
%                 
% %             figure(3)
% %             histogram(IDX,'normalization','probability')
% %                 
% %             figure(4)
% %             plot(time_all(1:time_set), IDX,'.-')  % plot(time_all, IDX,'.-')
% %                 
%             f1 = figure(100);
% %             [dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test(rel_phase_w_mean(1,:)', chan_coord(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%             [dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test(rel_phase_w_mean(1,:)', chan_coord_yx(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%             close(f1)
%                     
%                 
%                 %% get centroid_A and centroid_B
%                 
%                 centroid_K = nan(size(temp1,1),size(temp1,2),K);
%                 
%                 for j=1:K
%                     for i=1:length(topo_idx_x)
%                         centroid_K(topo_idx_x(i),topo_idx_y(i),j) = C(j,i); 
%                     end
%                 end
%                 
%                 figure1_b=figure;
%                 set (figure1_b, 'resize','off')
%                 set (figure1_b, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 40 35];  % [0 0 40 24];
%                 set(figure1_b, 'Position', pos, 'Units', 'centimeters');
%                 for j=1:K
%                     subplot(4,5,j)
%                     topoplot_figure(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
%                     noK=length(find(IDX==j));
%                     title(['Cluster ' num2str(j) ' : ' num2str(noK)])
%                 end
%                 print(figure1_b,'-dtiffn', '-r300', [outputfolder '\figures\Kmean cluster topo - ' new_name])
%                 close gcf 
% 
%                     %figure(1)
%                     %imagesc(centroid_A);set(gca,'YDir','normal');
%                     %colormap(jet);
%                     %axis square;
%                     
%                     %figure(2)
%                     %imagesc(centroid_B);set(gca,'YDir','normal');
%                     %colormap(jet);
%                     %axis square;
%                 
%                 %% Use TSNE algorithm to see if the clustering is well performed
%                 
% %                 tic
%                 tsne_topo = tsne(topo_vector);
% 
%                 figure1_c=figure;
%                 set (figure1_c, 'resize','off')
%                 set (figure1_c, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_c, 'Position', pos, 'Units', 'centimeters');%                 
%                 gscatter(tsne_topo(:,1),tsne_topo(:,2),IDX)
%                 print(figure1_c,'-dtiffn', '-r300', [outputfolder '\figures\Kmean tsne - ' new_name])
%                 close gcf 
% %                 toc
%                     
%                 % tic
%                 Y = pdist(topo_vector);
%                 YS = squareform(Y);
%                 Z = linkage(Y);
%                 figure1_d=figure;
%                 set (figure1_d, 'resize','off')
%                 set (figure1_d, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_d, 'Position', pos, 'Units', 'centimeters');%     
%                 dendrogram(Z);
%                 print(figure1_d,'-dtiffn', '-r300', [outputfolder '\figures\Dendrogram - ' new_name])
%                 close gcf 
%                 % toc
%                 
%                 %% dPLI and dwPLI computation, and ploting topoplot based on dPLI and dwPLI
%                 
%                 % whole time
%                 dPLI=d_PhaseLagIndex2(H);
%                 PLI = abs(dPLI);
%                 PLI(1:1+size(PLI,1):end) = nan;
%                 fC_PLI = nanmean(PLI);
%                 
%                 dwPLI=wd_PhaseLagIndex3_f(H); % dw_PhaseLagIndex(H);
%                 wPLI = abs(dwPLI);
%                 wPLI(1:1+size(wPLI,1):end) = nan;
%                 fC_wPLI = nanmean(wPLI);
%                 
%                 figure1_e=figure;
%                 set (figure1_e, 'resize','off')
%                 set (figure1_e, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_e, 'Position', pos, 'Units', 'centimeters');%   
% 
%                 topoplot_general_test(fC_PLI, chan_coord_yx(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 title('Functional Connectivity by PLI')
%                 print(figure1_e,'-dtiffn', '-r300', [outputfolder '\figures\FC PLI - ' new_name])
%                 close gcf 
% 
%                 figure1_f=figure;
%                 set (figure1_f, 'resize','off')
%                 set (figure1_f, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_f, 'Position', pos, 'Units', 'centimeters');%   
% 
%                 topoplot_general_test(fC_wPLI, chan_coord_yx(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 title('Functional Connectivity by wPLI')
%                 print(figure1_f,'-dtiffn', '-r300', [outputfolder '\figures\FC wPLI - ' new_name])
%                 close gcf 
% 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % 20 sec
%                 dPLI_20s=d_PhaseLagIndex2(H(end-sf*60+1:end-sf*40,:));
%                 PLI_20s = abs(dPLI);
%                 PLI_20s(1:1+size(PLI_20s,1):end) = nan;
%                 fC_PLI_20s = nanmean(PLI_20s);
%                 
%                 dwPLI_20s=wd_PhaseLagIndex3_f(H(end-sf*60+1:end-sf*40,:)); % dw_PhaseLagIndex(H);
%                 wPLI_20s = abs(dwPLI);
%                 wPLI_20s(1:1+size(wPLI_20s,1):end) = nan;
%                 fC_wPLI_20s = nanmean(wPLI_20s);
%                 
%                 figure1_e=figure;
%                 set (figure1_e, 'resize','off')
%                 set (figure1_e, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_e, 'Position', pos, 'Units', 'centimeters');%   
% 
%                 topoplot_general_test(fC_PLI_20s, chan_coord_yx(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 title('Functional Connectivity by PLI (20s)')
%                 print(figure1_e,'-dtiffn', '-r300', [outputfolder '\figures\FC PLI (20s) - ' new_name])
%                 close gcf 
% 
%                 figure1_f=figure;
%                 set (figure1_f, 'resize','off')
%                 set (figure1_f, 'Units', 'centimeters', 'Color', 'white');   % whole figure size control
%                 pos = [0 0 20 17];  % [0 0 40 24];
%                 set(figure1_f, 'Position', pos, 'Units', 'centimeters');%   
% 
%                 topoplot_general_test(fC_wPLI_20s, chan_coord_yx(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 title('Functional Connectivity by wPLI (20s)')
%                 print(figure1_f,'-dtiffn', '-r300', [outputfolder '\figures\FC wPLI (20s) - ' new_name])
%                 close gcf 
%                 
%                 %% Looking at the time series of "internal" vs "external" fluctuation 
%                 
%                 corr_pf = zeros(1,size(rel_phase_w_mean,1));
%                 for i=1:size(rel_phase_w_mean,1)
%                    corr_pf(i) = corr(rel_phase_w_mean(i,:)',fC_wPLI','type','Spearman');
%                 end
%                 
% %                 figure(8)
% %                 plot(time_all(1:length(corr_pf)),corr_pf)
% %                 xlabel('time(s)')
% %                 ylabel('corr(phase,fC)')
% %                 title('Time series for internal/external switching, by corr(phase,wPLI)')
% %                 
%                 corr_pf_sign = sign(corr_pf);
%                 
% %                 figure(9)
% %                 plot(time_all(1:length(corr_pf)),corr_pf_sign)
% %                 xlabel('time(s)')
% %                 ylabel('sign of corr(phase,fC)')
% %                 title('Binary time series for internal/external switching, by sign( corr(phase,wPLI) )')
%                 
% %                 idx_norm = (IDX.*2) -3; % (IDX-2.5)/1.5;
% %                 figure(10)
% %                 plot(time_all(1:length(idx_norm)),idx_norm, 'b.-'); hold on; 
% %                 plot(time_all(1:length(idx_norm)),corr_pf_sign(1:length(idx_norm)), 'r.-');
% %                 corr_kmean_pfsign = corr(corr_pf_sign(1:length(idx_norm))', idx_norm);
% %                 legend(['by K-mean clustering' ,'by sign( corr(phase,wPLI) )'])
% %                 title('Comparing two time series for internal/external switching')
% 
% 
%         save([ outputfolder new_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_alp_K_' num2str(K) ], ... 
%               'filename','filepath','pz_ox','poz_ox','plot_ch_idx','b','P','tr_mk','chan_coord_yx','rel_phase_w_mean','time_pt', ...
%               'topo','topo_vector','eva*','K','IDX','C','SUMD','D','sil','SIL','centroid_K','borderCoords','xx','yy','Coord', ...
%               'tsne_topo','Z','*PLI*','corr_pf','corr_pf_sign', '-v7.3');
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        display([num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm'])
        display(['Subject ' num2str(sub_no) '/22, ' 'session ' num2str(ses_no) '/' num2str(length(ses_list)) ', state ' num2str(st_no) '/' num2str(length(st_list))])
        real_time=clock;
        display([num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm'])

        end % st_no
        

    end % ses_no
    

end % sub_no


%%
% making topo_vector

clear all
close all
clc
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% loading file names in data directory
datapath=['D:\Slowave Dropbox\Kim Hyoungkyu\SKKU\NKI_download\raw_data\'];
datafold=dir([datapath 'sub-*']);
tic
real_time=clock;
[num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm']
bandname={'d', 't', 'a1', 'a2', 'b', 'r', 'w', '10'};
bandname2={'DELTA', 'THETA', 'ALPHA', 'ALPHA_narrow', 'BETA', 'GAMMA', 'WHOLE', '10Hz', 'ALPHA', 'Wave10'};
band_freq = [0.5 4 ; 4 8 ; 8 12 ; 9 11 ; 15 30 ; 30 50; 1 50; 9.5 10.5; 8 13];   % spectrogram paper 0.05~50Hz
electrodeExclude = 32; % ECG channel
%%%%%%%%%%%%%%%%%%%%
band_no=3;  
%%%%%%%%%%%%%%%%%%%%
sf=5000;  % sampling frequency = 500 Hz
% width=5;  % for wavelet filter parameter
% Nsurro=20;
% tr_drop=[0 1 2];  % 0 : no trial drop, 1 : light trial drop, 2 : heavy channel drop 
ch_no_F=sort([1 2 60 38 39 46 3 32 17 33 4 47]); % Fp1 Fp2 Fp AF3 AF4 F5 F3 F1 Fz F2 F4 F6 
ch_no_P=sort([50 7 36 19 37 8 51]);    % P5 P3 P1 Pz P2 P4 P6 
ch_no_O=sort([44 31 45]);              % PO3 POz Po4  //  9 10 20 - O1 O2 Oz 
ch_Cz=sort([34 18 35 23 61 24]);       % C1 Cz C2 CP1 CPz CP2

reference_name={'Single', 'Average', 'Laplacian', 'Mastoid'};  %, 'AverageTB', 'LaplacianTB'};
ref_type_no=2;  % 1,2,3,4
ANCFlag=1; % 1 or 0
outputfolder_eeg = ['D:\Slowave Dropbox\Kim Hyoungkyu\SKKU\NKI_download\preproc_self\eeg\'];

No_Pz = [];   pz_idx = [];   poz_idx = [];
for sub_no=1:22 % [18 8 19 9 20 10 21 11 22 12] % 1:length(datafold) %  1:7 8:12 13:17 18:22
    ses_list = dir([datapath datafold(sub_no).name '\ses-*']);

    for ses_no=1 % 1:length(ses_list)

        st_list = dir([datapath datafold(sub_no).name '\' ses_list(ses_no).name '\eeg\*.set']);

        for st_no=1:length(st_list)
            filename = st_list(st_no).name;
            filename_tmp = filename;
            filename_tmp = strrep(filename_tmp,'.set','');
            filename = [filename_tmp '_GAPArm.set'];
            
            filepath = ['D:\Slowave Dropbox\Kim Hyoungkyu\SKKU\NKI_download\preproc_self\eeg'];
            
            if exist([filepath '\' filename]) == 2 && contains(filename, 'rest')

            %%%%%%%%%%%%%%%%%%%%% name of the current file to use
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            EEG = pop_loadset('filename',filename,'filepath',filepath);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


            new_name = strrep([filename(1:end-4)],'_',' ');

            %%%%%%%%%%%%%%%%%%%%% FILTERING STEP
%                 band_no=4;
                    EEG = pop_firws(EEG, 'fcutoff', band_freq(band_no,1), 'ftype', 'highpass', 'wtype', 'blackman', 'forder', 2750, 'minphase', 0);
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
                    EEG = pop_firws(EEG, 'fcutoff', band_freq(band_no,2), 'ftype', 'lowpass', 'wtype', 'blackman', 'forder', 2750, 'minphase', 0);
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

            data=EEG.data'; % double(EEG.data)';

            tr_mk=[];
            mk_count=0;
            for mk_no=1:length(EEG.event)
                marker_tmp = EEG.event(mk_no).type;
                if  marker_tmp(1) == 'R'
                    mk_count = mk_count +1;
                    tr_mk(mk_count,1) = EEG.event(mk_no).latency;
                end % if
            end % mk_no

            if length(tr_mk)<1
                mk_count=0;
                for mk_no=1:length(EEG.event)
                    marker_tmp = EEG.event(mk_no).type;
                    if  marker_tmp(1) == 'q'
                        mk_count = mk_count +1;
                        tr_mk(mk_count,1) = EEG.event(mk_no).latency;
                    end % if
                end % mk_no
                tr_mk=tr_mk(1:end-4,1);
            end

            smooth = 20;
            tend = 1; % 1/3; % portion of the time which you want to make the movie with
                
            % Define file name and load the file.
            H = data(round(tr_mk(2)+1):round(tr_mk(end))+10500,:);
            H = downsample(H,20);
            for ch_no=1:length(EEG.chanlocs)
                chan_coord(ch_no,1)=EEG.chanlocs(ch_no).X;
                chan_coord(ch_no,2)=EEG.chanlocs(ch_no).Y;
                chan_coord(ch_no,3)=EEG.chanlocs(ch_no).Z;
            end % ch_no
            chan_coord_yx = [ chan_coord(:,2) chan_coord(:,1) chan_coord(:,3) ];
                
            % Compute relative phase
            [rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta] = cal_rel_phase_v2(H) ;
                
            % Moving time window
            time_moving =5; % 25; % 100ms for 250Hz SF
            time_window =5; % 25; % 1050; % 50;
            dt=1;
            Fs = 250; % 5000;
                
            [ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);
                
            rel_p =  double(rel_phase_w_mean);

%             load_name = dir([ outputfolder new_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_alp_K_*' ]);
%             load([ datapath datafold(sub_no).name '\' ses_list(ses_no).name '\' load_name.name ], 'tr_mk','chan_coord_yx','time_pt', 'topo','topo_vector','borderCoords','xx','yy','Coord');
%                 
            % making topo and movie
%             vidfile = VideoWriter([ outputfolder '\figures\' new_name '_tm' num2str(time_moving) '_tw' num2str(time_window) ... 
%                                   '_sm' num2str(smooth) '_' num2str(band_freq(band_no,1)) 'to' num2str(band_freq(band_no,2)) 'Hz'] ,'MPEG-4');
%             vidfile.FrameRate = round(Fs/time_moving/dt); % 0.48 for 2100ms
%             vidfile.Quality = 75;
%             open(vidfile);    
%             FrameRate = round(Fs/time_moving/dt);
%                 
            time_dt = 1/Fs*time_moving;
            time_now = 0 + (time_window)/2/Fs;
            time_all = [1:dt:round(size(rel_p,1)./tend )]';
                
                
            formatSpec = '%.2f';
            time_set = size(time_all,1); % 100;
            topo = cell(time_set,1);  % cell(size(time_all,1),1);
                
            tic
            warning('off','all')
            for i=1:dt:round(size(rel_p,1)/tend)
                    
%                 fig_H=figure;
                topo{i} = topoplot_general_test_without_figure(rel_p(i,:)', chan_coord_yx(:,1:2),'smooth',smooth,'scatter',1);
%                 title(  [ '\rm' num2str(time_now,formatSpec) ' S'; ] ,'fontsize',16 )
%                 drawnow
%                 time_all(i)=time_now;
%                 F(i) = getframe(gcf); 
%                 writeVideo(vidfile,F(i));
%                 time_now = time_now + dt*time_dt;
%                 close
            end
%             close(vidfile)
%             warning('on','all')
            toc


%                 %% dPLI and dwPLI computation
%                 
%                 % whole time
% 
%                 for tr_no = 1:(length(H)/250-10.5)/2.1
%                     data_tr = H((tr_no-1)*(250*2.1)+1:(tr_no-1)*(250*2.1)+(10.5*250) , :);
%                     dPLI(:,:,tr_no)=d_PhaseLagIndex2(data_tr);
%                     dwPLI(:,:,tr_no)=wd_PhaseLagIndex3_f(data_tr);
%                 end % tr_no


        save([filepath '\analysis\topo_20ms_wo_reref\' new_name ' 20ms_topo_wo_reref'], ... 
              'filename','filepath','tr_mk','rel_phase_w_mean','time_pt','topo', '-v7.3');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        display([num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm'])
        display(['Subject ' num2str(sub_no) '/22, ' 'session ' num2str(ses_no) '/' num2str(length(ses_list)) ', state ' num2str(st_no) '/' num2str(length(st_list))])
        real_time=clock;
        display([num2str(real_time(1)) '.' num2str(real_time(2)) '.' num2str(real_time(3)) ' : ' num2str(real_time(4)) 'h' num2str(real_time(5)) 'm'])

        end % if exist([filepath '\' filename])

        end % st_no
        

    end % ses_no
    

end % sub_no


%%
% 
% % K=2, 100ms 
% clear all
% file_path='E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\';
% file_rest500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*500_*alp_K_2.mat']);
% file_chkr500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*500_*alp_K_2.mat']);
% file_incp500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*500_*alp_K_2.mat']);
% 
% file_rest10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*10500_*alp.mat']);
% file_chkr10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*10500_*alp.mat']);
% file_incp10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*10500_*alp.mat']);
% 
% % Kmean (IDX) - 500 (100ms)
% Kmean_rest_2100ms=nan(288,3);
% Kmean_chkr_2100ms=nan(103,3);
% Kmean_incp_2100ms=nan(293,3);
% for session_no=1:3
%     load ([file_path file_rest10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_2100ms(:,session_no) = IDX;
% %     IDX_ratio(1,session_no) = find(IDX==1)/
%     load ([file_path file_chkr10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_2100ms(1:length(IDX),session_no) = IDX;
%     load ([file_path file_incp10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_2100ms(:,session_no) = IDX;
%     
% end % session_no
% 
% % Kmean (IDX) - 500 (100ms)
% Kmean_rest_100ms=nan(288*21,3);
% Kmean_chkr_100ms=nan(103*21,3);
% Kmean_incp_100ms=nan(293*21,3);
% 
% Kmean_rest_2100ms_from_100ms=nan(288,3);
% Kmean_chkr_2100ms_from_100ms=nan(103,3);
% Kmean_incp_2100ms_from_100ms=nan(293,3);
% 
% for session_no=1:3
%     load ([file_path file_rest500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_100ms(:,session_no) = IDX;
%     for i=1:length(IDX)/21
%         Kmean_rest_2100ms_from_100ms(i,session_no) = round(mean(IDX( (i-1)*21+1:i*21 )));
%     end 
%     [Kmean_corr_r(1,session_no), Kmean_corr_p(1,session_no)]=corr(Kmean_rest_2100ms_from_100ms(:,session_no), Kmean_rest_2100ms(:,session_no));
%     if Kmean_corr_r(1,session_no)>0
%         Kmean_ratio1(1,session_no) = length(find(Kmean_rest_2100ms(:,session_no)==1)) / length(find(Kmean_rest_2100ms(:,session_no)==2));
%         Kmean_ratio2(1,session_no) = length(find(Kmean_rest_2100ms_from_100ms(:,session_no)==1)) / length(find(Kmean_rest_2100ms_from_100ms(:,session_no)==2));
%     elseif Kmean_corr_r(1,session_no)<0
%         Kmean_ratio1(1,session_no) = length(find(Kmean_rest_2100ms(:,session_no)==1)) / length(find(Kmean_rest_2100ms(:,session_no)==2));
%         Kmean_ratio2(1,session_no) = length(find(Kmean_rest_2100ms_from_100ms(:,session_no)==2)) / length(find(Kmean_rest_2100ms_from_100ms(:,session_no)==1));
%     end % if 
%     if Kmean_ratio1(1,session_no)>1
%         Kmean_ratio1(1,session_no)=1/Kmean_ratio1(1,session_no);
%     end
%     if Kmean_ratio2(1,session_no)>1
%         Kmean_ratio2(1,session_no)=1/Kmean_ratio2(1,session_no);
%     end
% 
%     load ([file_path file_chkr500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_100ms(1:length(IDX),session_no) = IDX;
%     for i=1:length(IDX)/21
%         Kmean_chkr_2100ms_from_100ms(i,session_no) = round(mean(IDX( (i-1)*21+1:i*21 )));
%     end 
%     [Kmean_corr_r(2,session_no), Kmean_corr_p(2,session_no)]=corr(Kmean_chkr_2100ms_from_100ms(1:98,session_no), Kmean_chkr_2100ms(1:98,session_no));
%     if Kmean_corr_r(2,session_no)>0
%         Kmean_ratio1(2,session_no) = length(find(Kmean_chkr_2100ms(:,session_no)==1)) / length(find(Kmean_chkr_2100ms(:,session_no)==2));
%         Kmean_ratio2(2,session_no) = length(find(Kmean_chkr_2100ms_from_100ms(:,session_no)==1)) / length(find(Kmean_chkr_2100ms_from_100ms(:,session_no)==2));
%     elseif Kmean_corr_r(2,session_no)<0
%         Kmean_ratio1(2,session_no) = length(find(Kmean_chkr_2100ms(:,session_no)==1)) / length(find(Kmean_chkr_2100ms(:,session_no)==2));
%         Kmean_ratio2(2,session_no) = length(find(Kmean_chkr_2100ms_from_100ms(:,session_no)==2)) / length(find(Kmean_chkr_2100ms_from_100ms(:,session_no)==1));
%     end % if 
%     if Kmean_ratio1(2,session_no)>1
%         Kmean_ratio1(2,session_no)=1/Kmean_ratio1(2,session_no);
%     end
%     if Kmean_ratio2(2,session_no)>1
%         Kmean_ratio2(2,session_no)=1/Kmean_ratio2(2,session_no);
%     end
% 
%     load ([file_path file_incp500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_100ms(:,session_no) = IDX;
%     for i=1:length(IDX)/21
%         Kmean_incp_2100ms_from_100ms(i,session_no) = round(mean(IDX( (i-1)*21+1:i*21 )));
%     end 
%     [Kmean_corr_r(3,session_no), Kmean_corr_p(3,session_no)]=corr(Kmean_incp_2100ms_from_100ms(:,session_no), Kmean_incp_2100ms(:,session_no));
%     if Kmean_corr_r(3,session_no)>0
%         Kmean_ratio1(3,session_no) = length(find(Kmean_incp_2100ms(:,session_no)==1)) / length(find(Kmean_incp_2100ms(:,session_no)==2));
%         Kmean_ratio2(3,session_no) = length(find(Kmean_incp_2100ms_from_100ms(:,session_no)==1)) / length(find(Kmean_incp_2100ms_from_100ms(:,session_no)==2));
%     elseif Kmean_corr_r(3,session_no)<0
%         Kmean_ratio1(3,session_no) = length(find(Kmean_incp_2100ms(:,session_no)==1)) / length(find(Kmean_incp_2100ms(:,session_no)==2));
%         Kmean_ratio2(3,session_no) = length(find(Kmean_incp_2100ms_from_100ms(:,session_no)==2)) / length(find(Kmean_incp_2100ms_from_100ms(:,session_no)==1));
%     end % if 
%     if Kmean_ratio1(3,session_no)>1
%         Kmean_ratio1(3,session_no)=1/Kmean_ratio1(3,session_no);
%     end
%     if Kmean_ratio2(3,session_no)>1
%         Kmean_ratio2(3,session_no)=1/Kmean_ratio2(3,session_no);
%     end
% 
% end % session_no
% 
% save([file_path 'K2_mean_IDX_2100ms_100ms_rest_checker_inscape_session123'], 'Kmean*')
% 
% 
% 
% %%
% load ([file_path 'topo_parameters'])
% smooth=20;
% shiftpreset = 0;
% 
%                 for j=1:2
%                     figure(100+j-2)
%                     topoplot_figure_tmp(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
% %                     topoplot_figure(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
% %                     topoplot_general_test(centroid_K(:,:,j), chan_coord(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 end
% 
% %%
% session_no=1;
% % load ([file_path file_rest10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% % IDX_tmp=IDX;
% % load ([file_path file_rest500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% % IDX_0p1s=IDX;
% 
% % load ([file_path file_chkr10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% % IDX_tmp=IDX;
% % load ([file_path file_chkr500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% % IDX_0p1s=IDX;
% 
% load ([file_path file_incp10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% IDX_tmp=IDX;
% load ([file_path file_incp500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
% IDX_0p1s=IDX;
% 
% IDX_sum=[];
% for i=1:length(IDX_tmp)
% %     IDX_sum(i,1) = round(mean(IDX_0p21s( (i-1)*10+1:i*10 )));
%     IDX_sum(i,1) = round(mean(IDX_0p1s( (i-1)*21+1:i*21 )));
% end 
% [r,p]=corr(IDX_tmp, IDX_sum*(-1), 'type', 'spearman')
% %%
% figure; 
% plot(IDX_tmp, '-', 'markersize', 24)
% hold on
% plot(IDX_sum, '-', 'markersize', 24)
% 
% %%
% 
% figure;
% gscatter(tsne_topo(:,1),tsne_topo(:,2),IDX)
% 
% 
% %%
% figure; 
% plot(smooth(IDX_1050,20)); hold on
% plot([1:10:2780], IDX(11:end))
% xlim([1 500])
% 
% %%
% 
% % K=2
% clear all
% file_path='E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\';
% file_rest1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*1050_*alp.mat']);
% file_chkr1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*1050_*alp.mat']);
% file_incp1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*1050_*alp.mat']);
% 
% file_rest10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*10500_*alp.mat']);
% file_chkr10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*10500_*alp.mat']);
% file_incp10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*10500_*alp.mat']);
% 
% % Kmean (IDX) - 10500 (2100ms)
% Kmean_rest_2100ms=nan(288,3);
% Kmean_chkr_2100ms=nan(103,3);
% Kmean_incp_2100ms=nan(293,3);
% for session_no=2 % 1:3
%     load ([file_path file_rest10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_2100ms(:,session_no) = IDX;
% %     IDX_ratio(1,session_no) = find(IDX==1)/
%     load ([file_path file_chkr10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_2100ms(1:length(IDX),session_no) = IDX;
%     load ([file_path file_incp10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_2100ms(:,session_no) = IDX;
%     
% end % session_no
% 
% % Kmean (IDX) - 1050 (210ms)
% Kmean_rest_210ms=nan(2880,3);
% Kmean_chkr_210ms=nan(1030,3);
% Kmean_incp_210ms=nan(2930,3);
% for session_no=2 % 1:3
%     load ([file_path file_rest1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_210ms(:,session_no) = IDX;
% %     IDX_ratio(1,session_no) = find(IDX==1)/
%     load ([file_path file_chkr1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_210ms(1:length(IDX),session_no) = IDX;
%     load ([file_path file_incp1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_210ms(:,session_no) = IDX;
%     
% end % session_no
% 
% 
% %%
% load ([file_path 'topo_parameters'])
% smooth=20;
% shiftpreset = 0;
% 
%                 for j=1:2
%                     figure(100+j-2)
%                     topoplot_figure_tmp(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
% %                     topoplot_figure_tmp(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
% %                     topoplot_general_test(centroid_K(:,:,j), chan_coord(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
%                 end
% 
% %%
% for i=1:288
%     IDX_sum(i,1) = round(mean(IDX_0p21s( (i-1)*10+1:i*10 )));
% end 
% %%
% 
% figure;
% gscatter(tsne_topo(:,1),tsne_topo(:,2),IDX)
% 
% 
% %%
% figure; 
% plot(smooth(IDX_1050,20)); hold on
% plot([1:10:2780], IDX(11:end))
% xlim([1 500])
% 
% %%
% 
% % K=4
% clear all
% file_path='E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\';
% file_rest1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*1050_*_K_4.mat']);
% file_chkr1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*1050_*_K_4.mat']);
% file_incp1050 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*1050_*_K_4.mat']);
% 
% file_rest10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*rest*10500_*_K_4.mat']);
% file_chkr10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*checker*10500_*_K_4.mat']);
% file_incp10500 = dir(['E:\Dropbox\Cha-Moon\NKI EEG-fMRI data\sample_data\analysis_results\*inscape*10500_*_K_4.mat']);
% 
% % Kmean (IDX) - 10500 (2100ms)
% Kmean_rest_2100ms=nan(288,3);
% Kmean_chkr_2100ms=nan(103,3);
% Kmean_incp_2100ms=nan(293,3);
% for session_no=2 % 1:3
%     load ([file_path file_rest10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_2100ms(:,session_no) = IDX;
% %     IDX_ratio(1,session_no) = find(IDX==1)/
%     load ([file_path file_chkr10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_2100ms(1:length(IDX),session_no) = IDX;
%     load ([file_path file_incp10500(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_2100ms(:,session_no) = IDX;
%     
% end % session_no
% 
% % Kmean (IDX) - 1050 (210ms)
% Kmean_rest_210ms=nan(2880,3);
% Kmean_chkr_210ms=nan(1030,3);
% Kmean_incp_210ms=nan(2930,3);
% for session_no=2 % 1:3
%     load ([file_path file_rest1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_rest_210ms(:,session_no) = IDX;
% %     IDX_ratio(1,session_no) = find(IDX==1)/
%     load ([file_path file_chkr1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_chkr_210ms(1:length(IDX),session_no) = IDX;
%     load ([file_path file_incp1050(session_no).name], 'IDX', 'centroid_K', 'tsne_topo', 'eva')
%     Kmean_incp_210ms(:,session_no) = IDX;
%     
% end % session_no
% 
% 
% %%
% load ([file_path 'topo_parameters'])
%                 for j=1:4
%                     figure(100+j+30)
%                     topoplot_figure(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
%                 end
% 
% %%
% 
% figure;
% gscatter(tsne_topo(:,1),tsne_topo(:,2),IDX)
% 
% 
% %%
% figure; 
% plot(smooth(IDX_1050,20)); hold on
% plot([1:10:2780], IDX(11:end))
% xlim([1 500])
% 
% 









