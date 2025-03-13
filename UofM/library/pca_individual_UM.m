%% PCA 3-D, spatial PCA analysis 1
% clear;

% state = 'EC';
% filter = 'alpha';
time_moving = 10;
time_window = 10;
smooth = 20;

load_path = [pwd() '/..'];
save_path = [load_path '/results_yjp'];
load([save_path '/UM_info_AE.mat']); 
load('references/topo_vector_idx_v4.mat');

load(['../results_yjp/movie_rel_phase/band_alpha_re_whole/' subjname_ae{idx} '/' ...
    subjname_ae{idx} '_st_' state '_band_alpha_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo.mat'],'topo');
load(['../results_yjp/k_means/band_alpha_re_whole_each_mask/' subjname_ae{idx} '/' ...
    subjname_ae{idx} '_st_' state '_band_alpha_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans.mat'],'IDX');

topo_size = length(topo);
topo_vector = zeros(topo_size,length(topo_idx_x));

% topo_vector: time point x frame vector
for i=1:topo_size
    T = topo{i};
    topo_vector(i,:) = T(topo_vector_idx_v4);
end

% % to check the correlation k-mean algorithm
% topo_vector = topo_vector - mean(topo_vector,2);
% Xnorm = sqrt(sum(topo_vector.^2,2));
% topo_vector = topo_vector./Xnorm;

tic;
[coeff,score,latent,tsquared,explained,mu] = pca(topo_vector);
toc;

if matches(state, 'EO')
    score1 = score;
    explained1 = explained;
else
    EO = load(['../results_yjp/k_means/band_alpha_re_whole_each_mask/' subjname_ae{idx} '/' ...
        subjname_ae{idx} '_st_EO_band_alpha_tm' num2str(time_moving) ...
        '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans_pca.mat'], 'coeff');
    score1 = topo_vector*EO.coeff;
    explained1 = var(score1,[],1)';
    explained1 = explained1/sum(explained1)*100;
end

save(['../results_yjp/k_means/band_alpha_re_whole_each_mask/' subjname_ae{idx} '/' ...
    subjname_ae{idx} '_st_' state '_band_alpha_tm' num2str(time_moving) ...
    '_tw' num2str(time_window) '_sm' num2str(smooth) '_topo_kmeans_pca.mat'], ...
    'coeff', 'score', 'score1', 'latent', 'tsquared', 'explained', 'explained1', 'mu', 'IDX', '-v7.3');

%% PCA 3-D, spatial PCA analysis test 4

fig = figure(511); clf;
fig.Position = [50, 100, 2500, 600];
subplot(1,3,1); grid on; hold on;
plot(explained1,'o-','LineWidth',2,'MarkerSize',10);
xlabel('Principal compnent index');
ylabel('Percentage');
title({'Percentage of' 'the total variance explained'}, 'FontSize',20);
xlim([1 15]);
ylim([0 60]);
set(gca, 'FontSize', 15);

explained_cdf = zeros(1,length(explained1));
explained_cdf(1) = explained1(1);
for h=2:length(explained1)
    explained_cdf(h)=explained1(h)+explained_cdf(h-1);
end

subplot(1,3,2); grid on; hold on;
plot(explained_cdf,'o-','LineWidth',2,'MarkerSize',10)
xlabel('principal compnent index')
ylabel('cumulative percentage')
title('cumulative percentage of the total variance explained')
xlim([1 28])
ylim([0 100])


subplot(1,3,3); grid on; hold on;
for i=1:max(IDX)
    temp_idx = find(IDX==i);
    b(i) = plot3(score1(temp_idx,1), score1(temp_idx,2), score1(temp_idx,3),'o', ...
        'MarkerSize',3);
end

xlabel('PC1');ylabel('PC2');zlabel('PC3')
view(375, 20);
axis equal;
xL = [-180 180]; yL = [-150 150]; zL = [-100 100];
xlim(xL); ylim(yL); zlim(zL);
line([0 0], yL, 'color', 'k' , 'linewidth', 2);  %x-axis
line(xL, [0 0], 'color', 'k' , 'linewidth', 2); 
line([0 0 ], [0 0], zL, 'color', 'k' , 'linewidth', 2 );
legend(b, {'1','2','3','4'}, 'Location','northeast');
