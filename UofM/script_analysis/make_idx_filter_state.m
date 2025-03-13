%% Readme
% Commnad of making script files
% mrun "make"
% example, mrun "make"

%% make directories when they dosen't exist
if ~exist("scripts", "dir")
    mkdir("scripts");
end
if ~exist("err", "dir")
    mkdir("err");
end
if ~exist("out", "dir")
    mkdir("out");
end

%% generate both script files and run file
fp = fopen("run.sh", "w");
mrun = ['/usr/local/MATLAB/R2021b/bin/matlab ' ...
    '-nodesktop -nosplash -singleCompThread ' ...
    '-batch'];

title = 'UM';
file_path = [pwd() '/..'];
save_path = [file_path '/results_yjp'];

load([save_path '/UM_info_ch128.mat']); 
% load([save_path '/UM_info_AE.mat']); 
% n_sub    = n_sub_ae;
% subjname = subjname_ae;
% labels   = labels_ae;
% load([save_path '/UM_info_BS.mat']); 
% n_sub    = n_sub_bs;
% subjname = subjname_bs;
% labels   = labels_bs;

filters = filter_names;
% filters = filter_names(:,[4]);
state_list = {'EO', 'EC', 'LOC', 'AN', 'BS', 'DS', 'DA', 'ROC'};
for idx = 1:n_sub
%     for fi = 2:7
    for fi = 4 %5:7%4
%         filter = filters{fi};

        for st = 1:length(state_list)
            try
                state = state_list{st};
    
                fname = ['scripts/' subjname{idx} '_b_' num2str(fi) '_st_' num2str(st) '_1.sh'];
                fprintf(fp, 'qsub %s\n', fname);
            
                bashfp = fopen(fname, "w");
                jobname = [subjname{idx} '_b_' num2str(fi) '_st_' num2str(st) '_1'];
            
                fprintf(bashfp, '#PBS -S /bin/bash\n');
                fprintf(bashfp, '#PBS -N %s\n',         jobname); % the name of job
                fprintf(bashfp, '#PBS -o out/%s.out\n', jobname); % output file
                fprintf(bashfp, '#PBS -e err/%s.err\n', jobname); % error log file
%                 fprintf(bashfp, '#PBS -l nodes=1\n');
                fprintf(bashfp, 'cd $PBS_O_WORKDIR\n\n');
            
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'preprocessing_slicing_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'make_movie_v2_rel_phase_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'Copy_of_make_movie_v2_rel_phase_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'calculate_kmeans_topo_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'k_means_clustering_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'Copy_of_k_means_clustering_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'k_means_clustering_w_EO_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'Copy_of_k_means_clustering_w_EO_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'k_means_clustering_w_comb_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'k_means_clustering_w_comb_w_NOTA_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'k_means_clustering_w_EC_w_noto_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'measure_ApEn_dwelling_time_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'measure_ApEn_dwelling_time_v2_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'measure_ApEn_dwelling_time_v3_nota_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'pca_individual_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'Copy_of_pca_individual_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; filter=''%s''; state=''%s''; ' ...
%                     'draw_topomaps_UM"'], mrun, idx, filter, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; ' ...
%                     'pca_test_UM"'], mrun, idx, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; ' ...
%                     'make_movie_v2_rel_phase_v2_UM"'], mrun, idx, state);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; filter=%d; ' ...
%                     'make_movie_v2_rel_phase_v2_UM_tmp"'], mrun, idx, state, fi);
                fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; filter=%d; ' ...
                    'cal_regression_w_individual_UM"'], mrun, idx, state, fi);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; filter=%d; ' ...
%                     'cal_regression_w_individual_w_bandmask_UM"'], mrun, idx, state, fi);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; filter=%d; ' ...
%                     'cal_regression_idx_prop_w_individual_UM"'], mrun, idx, state, fi);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; filter=%d; ' ...
%                     'cal_regression_idx_prop_w_individual_w_1min_UM"'], mrun, idx, state, fi);
%                 fprintf(bashfp, ['%s "user_addpath(false,false); idx=%d; state=''%s''; ' ...
%                     'pca_individual_UM_240619"'], mrun, idx, state);
                fclose(bashfp);
            end % try
        end % state
    end % filters
end % idx
fclose(fp);
