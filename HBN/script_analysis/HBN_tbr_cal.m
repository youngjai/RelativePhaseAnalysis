%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load EEG files
% 2. Extract data and event
% 3. Calculate theta (4–8 Hz) and beta (12–20 Hz) power for each segment
% 4. Compute TBR (theta/beta ratio) 
% 5. Plot figures

% 2024. 06. 10. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TBR Calculation
% Beta band: 12–20 Hz

file_path = 'D:\HealthyBrainNetwork_above11\EEGDownload4\C_w_124_bp_bc_original_new';
file_list = dir(fullfile(file_path, '*.mat'));

load('D:\HealthyBrainNetwork_above11\EEGDownload4\bc_idx_C_original_new_25p_last_93.mat');

Fs = 500;
Fz = 'E11';
Pz = 'E62';

for subj_no = 1:length(file_list)
    cd(file_path)

    input_fileName = file_list(subj_no).name;
    fileName = erase(input_fileName, '.mat');

    load(input_fileName);

    % extract EEG data and sampling rate
    data = EEG.data;
    Fs = EEG.srate;

    % extract event information
    event = EEG.event;
    event_sample = [event.sample];
    event_type = transpose({event(:).type});
    event_time = event_sample * Fs / 60;

    RS_total.event(subj_no).type = event_type;
    RS_total.event(subj_no).sample = event_sample;
    RS_total.event(subj_no).time = event_time;

    % initialize data structures
    RS_total.data(subj_no).epoch = {};
    RS_total.data(subj_no).epoch_fix = {};
    RS_total.data(subj_no).EO = {};
    RS_total.data(subj_no).EC = {};
    RS_total.state(subj_no).EO = [];
    RS_total.state(subj_no).EC = [];
    RS_total.state(subj_no).EO_event = [];
    RS_total.state(subj_no).EC_event = [];

    % align event samples
    for es = 1:11
        event_sample_new(es) = event_sample(es) - (event_sample(1) - 1);
    end
    event_sample = event_sample_new';

    % epoching and state segmentation
    for ch = 1:size(data, 1)
        for es = 1:10
            temp = nonzeros(data(ch, event_sample(es):event_sample(es + 1)));
            RS_total.data(subj_no).epoch(ch).dat(es, 1:length(temp)) = temp;
        end

        for es = 1:10
            temp = RS_total.data(subj_no).epoch(ch).dat(es, 1:20001);
            RS_total.data(subj_no).epoch_fix(ch).dat(es, 1:length(temp)) = temp;
        end

        for en = 1:10
            if mod(en, 2) == 0  % ec
                idx = en / 2;
                RS_total.data(subj_no).EC(ch).dat(idx, :) = ...
                    RS_total.data(subj_no).epoch_fix(ch).dat(en, 1:20001);
                temp2 = transpose(nonzeros(RS_total.data(subj_no).EC(ch).dat(idx, :)));

                if idx == 1
                    RS_total.state(subj_no).EC(ch, 1:length(temp2)) = temp2;
                    RS_total.state(subj_no).EC_event(idx) = length(temp2);
                else
                    temp3 = length(nonzeros(RS_total.state(subj_no).EC(ch, :)));
                    RS_total.state(subj_no).EC(ch, temp3+1:temp3+length(temp2)) = temp2;
                    RS_total.state(subj_no).EC_event(idx) = length(RS_total.state(subj_no).EC(ch, :));
                end
            else  % eo
                idx = (en + 1) / 2;
                RS_total.data(subj_no).EO(ch).dat(idx, :) = ...
                    RS_total.data(subj_no).epoch_fix(ch).dat(en, 1:10001);
                temp2 = transpose(nonzeros(RS_total.data(subj_no).EO(ch).dat(idx, :)));

                if idx == 1
                    RS_total.state(subj_no).EO(ch, 1:length(temp2)) = temp2;
                    RS_total.state(subj_no).EO_event(idx) = length(temp2);
                else
                    temp3 = length(nonzeros(RS_total.state(subj_no).EO(ch, :)));
                    RS_total.state(subj_no).EO(ch, temp3+1:temp3+length(temp2)) = temp2;
                    RS_total.state(subj_no).EO_event(idx) = length(RS_total.state(subj_no).EO(ch, :));
                end
            end
        end
    end

    % extract ec/eo data (clipped segment)
    ch_num = length(RS_total.data(subj_no).EO);
    for i = 1:5
        for j = 1:ch_num
            EC(i).dat(j, :) = RS_total.data(subj_no).EC(j).dat(i, 2701:19500);
            EO(i).dat(j, :) = RS_total.data(subj_no).EO(j).dat(i, 1101:9500);
        end
    end

    % TBR Computation
    TBR(subj_no).name = input_fileName(1:12);
    TBR(subj_no).ch_loc = bc_idx_last(subj_no).good;

    for j = 1:5
        eo_data = EO(j).dat;
        ec_data = EC(j).dat;

        TBR(subj_no).eo_theta{j} = bandpower(eo_data', Fs, [4 8]);
        TBR(subj_no).ec_theta{j} = bandpower(ec_data', Fs, [4 8]);
        eo_theta_all(j, :) = TBR(subj_no).eo_theta{j};
        ec_theta_all(j, :) = TBR(subj_no).ec_theta{j};

        TBR(subj_no).eo_beta{j} = bandpower(eo_data', Fs, [12 20]);
        TBR(subj_no).ec_beta{j} = bandpower(ec_data', Fs, [12 20]);
        eo_beta_all(j, :) = TBR(subj_no).eo_beta{j};
        ec_beta_all(j, :) = TBR(subj_no).ec_beta{j};
    end

    % average TBR across epochs
    TBR(subj_no).tbr_eo = mean(eo_theta_all) ./ mean(eo_beta_all);
    TBR(subj_no).tbr_ec = mean(ec_theta_all) ./ mean(ec_beta_all);

    % extract Fz and Pz channel indices
    sub_ch_loc = TBR(subj_no).ch_loc;
    for ch = 1:length(sub_ch_loc)
        label = sub_ch_loc(ch).labels;
        if strcmp(label, Fz)
            ch_num_Fz = ch;
        elseif strcmp(label, Pz)
            ch_num_Pz = ch;
        end
    end

    
    TBR(subj_no).Fz_eo = TBR(subj_no).tbr_eo(ch_num_Fz);
    TBR(subj_no).Fz_ec = TBR(subj_no).tbr_ec(ch_num_Fz);
    TBR(subj_no).Pz_eo = TBR(subj_no).tbr_eo(ch_num_Pz);
    TBR(subj_no).Pz_ec = TBR(subj_no).tbr_ec(ch_num_Pz);

    clear('EC', 'EO', 'eo*', 'ec*')
end

%% save Results
TBR_Control = TBR;
%TBR_Ain = TBR;

save('TBR_Control', 'TBR_Control');
% save('TBR_Ain', 'TBR_Ain');


%% figure 
TBR_Ain_Fz_ec_m = mean(TBR_Ain_Fz_ec);
TBR_Ain_Fz_eo_m = mean(TBR_Ain_Fz_eo);
TBR_Ain_Pz_ec_m = mean(TBR_Ain_Pz_ec);
TBR_Ain_Pz_eo_m = mean(TBR_Ain_Pz_eo);

TBR_C_Fz_ec_m = mean(TBR_C_Fz_ec);
TBR_C_Fz_eo_m = mean(TBR_C_Fz_eo);
TBR_C_Pz_ec_m = mean(TBR_C_Pz_ec);
TBR_C_Pz_eo_m = mean(TBR_C_Pz_eo);

TBR_Ain_Fz_ec_se = std(TBR_Ain_Fz_ec)/sqrt(length(TBR_Ain_Fz_ec));
TBR_Ain_Fz_eo_se = std(TBR_Ain_Fz_eo)/sqrt(length(TBR_Ain_Fz_eo));
TBR_Ain_Pz_ec_se = std(TBR_Ain_Pz_ec)/sqrt(length(TBR_Ain_Pz_ec));
TBR_Ain_Pz_eo_se = std(TBR_Ain_Pz_eo)/sqrt(length(TBR_Ain_Pz_eo));

TBR_C_Fz_ec_se = std(TBR_C_Fz_ec)/sqrt(length(TBR_C_Fz_ec));
TBR_C_Fz_eo_se = std(TBR_C_Fz_eo)/sqrt(length(TBR_C_Fz_eo));
TBR_C_Pz_ec_se = std(TBR_C_Pz_ec)/sqrt(length(TBR_C_Pz_ec));
TBR_C_Pz_eo_se = std(TBR_C_Pz_eo)/sqrt(length(TBR_C_Pz_eo));

Fz_Pz_se = horzcat(TBR_C_Fz_eo_se, TBR_Ain_Fz_eo_se, TBR_C_Fz_ec_se, TBR_Ain_Fz_ec_se, TBR_C_Pz_eo_se, TBR_Ain_Pz_eo_se, TBR_C_Pz_ec_se, TBR_Ain_Pz_ec_se);

Fz_eo = horzcat(TBR_C_Fz_eo_m, TBR_Ain_Fz_eo_m);
Fz_ec = horzcat(TBR_C_Fz_ec_m, TBR_Ain_Fz_ec_m);

Pz_eo = horzcat(TBR_C_Pz_eo_m, TBR_Ain_Pz_eo_m);
Pz_ec = horzcat(TBR_C_Pz_ec_m, TBR_Ain_Pz_ec_m);

Fz_Pz = horzcat(Fz_eo, Fz_ec, Pz_eo, Pz_ec);


[h_ec_fz, p_ec_fz] = ttest2(TBR_C_Fz_ec, TBR_Ain_Fz_ec);
[h_eo_fz, p_eo_fz] = ttest2(TBR_C_Fz_eo, TBR_Ain_Fz_eo);

[h_ec_pz, p_ec_pz] = ttest2(TBR_C_Pz_ec, TBR_Ain_Pz_ec);
[h_eo_pz, p_eo_pz] = ttest2(TBR_C_Pz_eo, TBR_Ain_Pz_eo);


fig = figure(1);
fig.Position = [50 50 800 500];

grid on
hold on
b1 = bar([1:2], Fz_eo);
b2 = bar([4:5], Fz_ec);
b3 = bar([8:9], Pz_eo);
b4 = bar([11:12], Pz_ec);

set(b1,'FaceColor', '#26A69A', 'FaceAlpha', 0.5);
set(b2,'FaceColor', '#26A69A', 'FaceAlpha', 0.6);
set(b3,'FaceColor', '#FF5722', 'FaceAlpha', 0.5);
set(b4,'FaceColor', '#FF5722', 'FaceAlpha', 0.6);

errorbar([1:2, 4:5, 8:9, 11:12], Fz_Pz, Fz_Pz_se, 'k.', 'LineWidth', 0.5);
set(gca, 'XTick', [1:2, 4:5, 8:9, 11:12], 'XTickLabel', {'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD', 'Control', 'ADHD'}, 'FontSize', 24);

ylabel('Theta/Beta Ratio', 'FontSize', 24);

% Shading for the intervals
shading_intervals = {[3, 6], [10, 13]};
y_limits = get(gca, 'YLim');  % Get the y-axis limits


% RED [0.984, 0.914, 0.906]

for i = 1:length(shading_intervals)
    x = shading_intervals{i};
    x_patch = [x(1), x(2), x(2), x(1)];  % x-coordinates for the shading rectangle
    y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];  % y-coordinates for the shading rectangle
    patch(x_patch, y_patch, [1, 0.95, 0.46], 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Add shading with light gray color
end

legend([b1, b3], {'Fz', 'Pz'}, 'FontSize', 12, 'location', 'northoutside');