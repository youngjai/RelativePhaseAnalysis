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
mrun = '/usr/local/MATLAB/R2021b/bin/matlab -nodesktop -nosplash -singleCompThread -batch';

title = 'UM';
file_path = [pwd() '/..'];
save_path = [file_path '/results_yjp'];

load([save_path '/UM_info_AE.mat']); 
n_sub = n_sub_ae;
subjname = subjname_ae;
% load([save_path '/UM_info_BS.mat']); 
% n_sub = n_sub_bs;
% subjname = subjname_bs;

filters = filter_names;
% measures = {'ApEn', 'Shannon', 'corrDim', 'std', 'PLE', 'imc_relp', 'MST'};
measures = {'PLE', 'ApEn', 'Shannon', 'corrDim', 'std', 'imc_relp'};
% measures = {'MST'};


for i = 1:length(filters)
    filter = filters{i};
    for j = 1:length(measures)
        measure = measures{j};

        if ~matches(measure, 'Shannon')
            % ---------- 5
            fname = ['scripts/' title '_idx_filter_' num2str(i) '_' num2str(j) '.sh'];
            fprintf(fp, 'qsub %s\n', fname);
            
            bashfp = fopen(fname, "w");
            jobname = [title '_' num2str(i) '_' num2str(j) '_d_' measure];
            
%             fprintf(bashfp, '#!/bin/bash\n');
%             fprintf(bashfp, '#SBATCH -J %s\n',         jobname); % the name of job
%             fprintf(bashfp, '#SBATCH -o out/%s.out\n', jobname); % output file
%             fprintf(bashfp, '#SBATCH -e err/%s.err\n', jobname); % error log file
%             fprintf(bashfp, '#SBATCH -t unlimited\n');
%             fprintf(bashfp, '#SBATCH -n 1\n');
            fprintf(bashfp, '#PBS -S /bin/bash\n');
            fprintf(bashfp, '#PBS -N %s\n',         jobname); % the name of job
            fprintf(bashfp, '#PBS -o out/%s.out\n', jobname); % output file
            fprintf(bashfp, '#PBS -e err/%s.err\n', jobname); % error log file
            fprintf(bashfp, '#PBS -l nodes=1\n');
            fprintf(bashfp, 'cd $PBS_O_WORKDIR\n\n');

            fprintf(bashfp, '%s "filter=''%s''; measure=''%s''; draw_time_scale_UM"', ...
                mrun, filter, measure);
            fclose(bashfp);
        else
            n_list = [2 5 10 20 50 100];
            for k = 1:length(n_list)
                n_bins = n_list(k);
                % ---------- 6
                fname = ['scripts/' title '_idx_filter_' num2str(i) '_' num2str(j) '_' num2str(k) '.sh'];
                fprintf(fp, 'qsub %s\n', fname);
            
                bashfp = fopen(fname, "w");
                jobname = [title '_' num2str(i) '_' num2str(j) '_' num2str(k) '_d_' measure];
            
%                 fprintf(bashfp, '#!/bin/bash\n');
%                 fprintf(bashfp, '#SBATCH -J %s\n',         jobname); % the name of job
%                 fprintf(bashfp, '#SBATCH -o out/%s.out\n', jobname); % output file
%                 fprintf(bashfp, '#SBATCH -e err/%s.err\n', jobname); % error log file
%                 fprintf(bashfp, '#SBATCH -t unlimited\n');
%                 fprintf(bashfp, '#SBATCH -n 1\n');
                fprintf(bashfp, '#PBS -S /bin/bash\n');
                fprintf(bashfp, '#PBS -N %s\n',         jobname); % the name of job
                fprintf(bashfp, '#PBS -o out/%s.out\n', jobname); % output file
                fprintf(bashfp, '#PBS -e err/%s.err\n', jobname); % error log file
                fprintf(bashfp, '#PBS -l nodes=1\n');
                fprintf(bashfp, 'cd $PBS_O_WORKDIR\n\n'); 

                fprintf(bashfp, '%s "filter=''%s''; measure=''Shannon''; n_bins=%d; draw_time_scale_UM"', ...
                    mrun, filter, n_bins);
                fclose(bashfp);
            end % n_bins
        end % Shannon or not
    end % measures
end % filters

