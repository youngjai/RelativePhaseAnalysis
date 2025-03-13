function [IDX,C,SUMD,D] = kmeans_with_mask_v2(mask, data, varargin)
% made by Youngjai Park, June 27th, 2023.
% To calculate the time series with considering a given mask
%
% mask : Canonical modes information, the size is K by topovector size.
% data : which is that we want to categorize as mask, the size is time by topovector size. 
% method : 'euclidean' or 'pearson', when we calculate the distance from mask, we use the method.

% set default values
method = 'euclidean';

try
    method = varargin{1};
end

data_size = size(data,1);
IDX = zeros([data_size,1]);
C = mask;
K = size(mask,1);
D = zeros([data_size,K]);

switch method
    case 'euclidean'
        for t = 1:data_size
            D(t,:) = vecnorm(C-data(t,:),2,2);
            [~,IDX(t)] = min(D(t,:));
        end
    case 'pearson'
        for t = 1:data_size
            for k = 1:K
                D(t,k) = corr(C(k,:)',data(t,:)');
            end
            [~,IDX(t)] = max(D(t,:));
        end
end
SUMD = sum(D,1)';

end