% put reference mask, it returns the new order with matching the mask
% for example: if you want to put mode #2 in the first place and mode #1 in
% the second place, new order=[2 1];
% K : k-mean number
% topo_size : put the topomap (i.e. topo{i})
% topo_vector_idx : e.g.) 31759 by 1 vector

% changed by YJ Park at Oct. 20th, 2023.

function [IDX, C, SUMD, D, rho] = change_cluster_idx_v3(K, topo_size, topo_vector_idx, ...
    IDX_old, C_old, SUMD_old, D_old)

    % make the mask
    size_x2 = round(topo_size/2); size_y2 = round(topo_size/2);
    mask_topo = nan(topo_size,topo_size,4);
    
    % front-to-back
    mask_topo(:,:,1) = -1; mask_topo(size_x2:end,:,1) = 1;
    % back-to-fornt
    mask_topo(:,:,4) = -1.*mask_topo(:,:,1);
    % left-to-right
    mask_topo(:,:,2) = -1; mask_topo(:,1:size_y2,2) = 1;
    % right-to-left
    mask_topo(:,:,3) = -1.*mask_topo(:,:,2);
    
    mask = zeros(K,length(topo_vector_idx));
    switch K
        case 2
            % front-to-back
            T = mask_topo(:,:,1); mask(1,:) = T(topo_vector_idx);
            % back-to-front
            T = mask_topo(:,:,4); mask(2,:) = T(topo_vector_idx);
        case 4
            for k = 1:K
                T = mask_topo(:,:,k); mask(k,:) = T(topo_vector_idx);
            end
    end


    % change the indexes
    rho = corr(mask', C_old');
    new_order = zeros([1 K]);
    for k = 1:K
        [~,new_idx] = max(rho(k,:));
        new_order(k) = new_idx;
    end

    old_order = 1:K;
    C(old_order,:) = C_old(new_order,:);
    SUMD(old_order,:) = SUMD_old(new_order,:);
    D(:,old_order) = D_old(:,new_order);
    rho(:,old_order) = rho(:,new_order);

    IDX = IDX_old;
    for i=1:K
        cio = old_order(i);
        cin = new_order(i);
        IDX(IDX_old==cin)=cio;
    end

end