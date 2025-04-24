
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load the combined cetroids
% 2. Assign directions to masks
% 3. Load all topo vector
% 4. Perform beta regression
% 5. Calculation properties
%
% 2024. 04. 09. Younghwa Cha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% beta regression

clear, close all;
mask = load('D:\Dropbox\Moonbrainlab\HBN\centroid_topo_ref\comb_centroids_20240311/combined_centroids_20240311.mat', 'C_comb').C_comb.EOEC;
load('D:\Dropbox\Moonbrainlab\HBN\centroid_topo_ref\comb_centroids_20240311/combined_centroids_20240311.mat', 'topo_idx*');
% x1 : centroid of C1 and -C4 (front lead)
% x2 : centroid of C2 and -C3 (left lead)
mask([3 4],:,:) = -1.*mask([3 4],:,:);
x1 = squeeze(mean(mask([1 4],:,:),1));
x1 = x1(topo_idx_comb); x1 = x1./norm(x1);
x2 = squeeze(mean(mask([2 3],:,:),1)); 
x2 = x2(topo_idx_comb); x2 = x2./norm(x2);
X = cat(2,x1,x2);


load('D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\topo_vector\Ain_eo_original_new_all_topo', 'all_topo_vector');
topo_vector = all_topo_vector(:,topo_idx_HBN2comb);
[n_T, n_pnt] = size(topo_vector);
beta = zeros([3,n_T]); % beta1, beta2, and epsilon
readme = ['x1 : centroid of C1 and -C4 (front lead); ' ...
    'x2 : centroid of C2 and -C3 (left lead); ' ...
    'beta: beta1, beta2, and epsilon;'];
y = topo_vector';

tic;
for t = 1:size(y,2)
    [beta(1:2,t),~,E] = mvregress(X,y(:,t));
    beta(3,t) = mean(E);
end
toc;

name = Ain_eo
savepath = ['D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result\' name];
save(savepath, 'beta', 'X', 'readme');


%% Calculation properties by subjects
% subject separtion
time = 336 %ec 336 / eo 168 (*100 ms) the number of frame
num = length(beta)/time/5;


for i = 1:num
    
    beta_sub{i} = beta(:, (time*5*(i-1)+1):(time*5*i));

end


% clustering 

K = 4; unit = 0.1; n_nulls = 10000;

for i = 1:num
    tic
    [IDX, w_prop] = cal_regression_clustering(beta_sub{i}, K, unit, n_nulls);
    toc
    
    IDX_sub{i} = IDX;
    w_prop_sub(i) = w_prop;

    prop = cal_transition_prop_v2(IDX, K, unit);
    prop_sub(i)= prop;
end

name = Ain_eo
savepath_cal_sub = ['D:\Dropbox\MoonBrainLab Raw Data\Healthy Brain Network data\regress_result\' name '_cal_sub'];
save(savepath_cal_sub, '*_sub');


