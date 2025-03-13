% made by Youngjai Park, 2023.03.29.

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];

load([save_path '/UM_info_BS.mat']);

fig = figure('NumberTitle','off', 'Visible','off');    % make frame without display
fig.Position(3:4) = [800,600];   
hold on;
for idx = 1:n_sub_bs
    subplot(6,1,idx);
    X = (timestemp_bs{2,idx}(:,2)-timestemp_bs{2,idx}(:,1))/fs;
    edges = 1:60;
    h = histogram(X,edges);
    x = flip(h.BinEdges(1:end-1));
    y = flip(h.Values);
    tot_time = 0;
    for i = 1:h.NumBins-1
        tot_time = tot_time + x(i)*y(i);
        if tot_time > 300
            break;
        end
    end
    tot_time = 0;
    y_cum = cumsum(y);
    for j = 1:h.NumBins-1
        tot_time = x(j)*y_cum(j);
        if tot_time > 300
            break;
        end
    end
    text(max(xlim)*0.9, max(ylim)*0.7, ...
        ['# of events: ' num2str(size(X,1))], ...
        'Horiz','center', 'Vert','middle');
    xline(x(i), LineWidth=5, Color='r', Label=['a) ' num2str(x(i))]);
    xline(x(j), LineWidth=5, Color='b', Label=['b) ' num2str(x(j))]);
    ylabel(subjname_bs{idx}, FontSize=12, Interpreter='none');
end
xlabel('Time (sec)', 'FontSize',14);
sgtitle('Histogram of burst suppression', 'FontSize',16);
exportgraphics(fig, [save_path '/burst_suppression/BS_histogram.png'], 'Resolution',600);
close(fig);
