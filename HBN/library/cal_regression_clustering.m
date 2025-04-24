function [IDX, w_prop] = cal_regression_clustering(beta, K, unit, n_nulls, IDX_in)
    if K == 4
        % indexing by using beta of regression
        X = beta(1,:)'; Y = beta(2,:)';
        IDX = 1.*((Y <X) & (Y>=-X)) ...
            + 2.*((Y>=X) & (Y>=-X)) ...
            + 3.*((Y <X) & (Y <-X)) ...
            + 4.*((Y>=X) & (Y <-X));
        weight = [X Y -Y -X];
    else
        X = beta(1,:)'; Y = beta(2,:)';
        IDX = IDX_in;
        weight = [X Y -Y -X sqrt(X.^2+Y.^2)];
    end


    w_prop = struct;
    w_prop.avg_weight = zeros([1,K]);
    w_prop.sum_weight = zeros([1,K]);
    w_prop.dwell_time = zeros([1,K]);
    w_prop.dwell_dist = cell([1,K]);
    w_prop.trans_freq_mat = zeros([K,K]);
    w_prop.trans_prob_mat = zeros([K,K]);
    w_prop.trans_zscore_mat = zeros([K,K]);
    w_prop.K = K;
    w_prop.unit = unit;
    w_prop.n_nulls = n_nulls;

    for k = 1:K
        IDX_1 = 1.*(IDX==k); % binarize the given time series
        % So, IDX_1 array has occuring transition event when an element of IDX_1 is 1 or -1.
             
        % calculate the mean weight of 'mode k'. (ratio)
        w_prop.avg_weight(k) = sum(weight(:,k).*(IDX_1))./length(IDX);

        % calculate the mean weight of each axis. (ratio)
        w_prop.sum_weight(k) = sum(weight(weight(:,k)>0,k))./length(IDX);
    
        % calculate the dwelling time. (Second)
        IDX_2 = diff([0; IDX_1; 0]); % add '0' at the first and the last
        trans_id = find(IDX_2~=0); % find indexes when transition occurs.
        w_prop.dwell_dist{k} = zeros([0.5*length(trans_id),1]);
        for t = 1:2:length(trans_id)
            w_prop.dwell_dist{k}(0.5*(t+1)) = ...
                sum(weight(trans_id(t):trans_id(t+1)-1,k))*unit;
        end
        w_prop.dwell_time(k) = mean(w_prop.dwell_dist{k});
    end

    % calculate a transition matrix per a second. (Hertz)
    IDX_trans = IDX(1:end-1) - IDX(2:end).^K;
    trans_wdist = vecnorm([X(1:end-1) Y(1:end-1)] ...
        - [X(2:end) Y(2:end)],2,2);
    trans_mat = zeros(K);
    for i = 1:K
        for j = 1:K
            % to make all cases as different values, we use IDX_trans
            % (i,j) means the transition from mode j to mode i 
            trans_mat(i,j) = sum(trans_wdist.*(IDX_trans==(j-i^K)));
        end
    end
    w_prop.trans_freq_mat = trans_mat./(unit*(length(IDX_trans)));
    w_prop.trans_prob_mat = trans_mat./sum(trans_mat,1);

    % calculate the transition matrix of null models by using random rotationg beta
    rnd = -pi + 2*pi*rand(n_nulls,1);
    trans_mat_null = zeros([K,K,n_nulls]);
    for r = 1:n_nulls
        rot = [cos(rnd(r)) -sin(rnd(r));
               sin(rnd(r))  cos(rnd(r))];
        beta_tmp = rot*beta(1:2,:);
        X = beta_tmp(1,:)'; Y = beta_tmp(2,:)';
        IDX_null = 1.*((Y <X) & (Y>=-X)) ...
            + 2.*((Y>=X) & (Y>=-X)) ...
            + 3.*((Y <X) & (Y <-X)) ...
            + 4.*((Y>=X) & (Y <-X));

        IDX_trans = IDX_null(1:end-1) - IDX_null(2:end).^K;
        trans_wdist = vecnorm([X(1:end-1) Y(1:end-1)] ...
            - [X(2:end) Y(2:end)],2,2);
        for i = 1:K
            for j = 1:K
                % to make all cases as different values, we use IDX_trans
                % (i,j) means the transition from mode j to mode i 
                trans_mat_null(i,j,r) = sum(trans_wdist.*(IDX_trans==(j-i^K)));
            end
        end
    end
    % z score by using null models
    w_prop.trans_zscore_mat = (trans_mat-mean(trans_mat_null,3))./std(trans_mat_null,0,3);

end

