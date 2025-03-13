function [beta, epsilon, readme] = cal_regression(mask, topo_vector)
    % x1 : centroid of C1 and -C4 (front lead)
    % x2 : centroid of C2 and -C3 (left lead)
    mask([3 4],:) = -1.*mask([3 4],:);
    x1 = mean(mask([1 4],:),1)'; x1 = x1./norm(x1);
    x2 = mean(mask([2 3],:),1)'; x2 = x2./norm(x2);
    X = cat(2,x1,x2);
    
    [n_T, n_pnt] = size(topo_vector);
    
    beta = zeros([2,n_T]); % beta1 and beta2
    epsilon = zeros([n_pnt,n_T]); % epsilon
    readme = ['x1 : centroid of C1 and -C4 (front lead); ' ...
        'x2 : centroid of C2 and -C3 (left lead); ' ...
        'beta: beta1 and beta2; ' ...
        'epsilon: epsilon;'];
    
    y = topo_vector';
    
    for t = 1:size(y,2)
        [beta(1:2,t),~,epsilon(:,t)] = mvregress(X,y(:,t));
    end

end