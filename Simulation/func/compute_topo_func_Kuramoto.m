
function [topo_out] = compute_topo_func_Kuramoto(X_prod, chan_coord_xy, param, SEED)
    if ~isfield(param,'last_50_flag'); param.last_50_flag = 0; end
    TSTART = tic;
    %% set params
    
    smooth = 80; % spatial smoothness
    tend_temp = 1; % portion of the time which you want to make the movie with
    tend = 1./tend_temp;

    % Moving time window
    time_moving = param.time_moving;  %  number of time samples, time_moving / Fs = seconds
    time_window = param.time_window;  %  number of time samples, time_window / Fs = seconds
    dt = 1;
%     Fs = 1000;         % frequency(Hz)

    %% Compute relative phase
    simulated_angle = X_prod(:,1:end-2);
    H = real(1.*exp(1i*simulated_angle));

    rel_phase_w = cal_rel_phase_v2(H) ;
    rel_phase_w_mean = moving_time_window(rel_phase_w, time_moving, time_window);
    rel_p = double(rel_phase_w_mean(1:round(size(rel_phase_w_mean,1)./tend),:));
    if param.last_50_flag
        numT = size(rel_p,1);
        % rel_p = rel_p(end-round(numT/2)+1:end,:);
        rel_p = rel_p(end-round(numT/3)*2+1:end,:);
    end
    time_all = 1:dt:round(size(rel_p,1)./tend);
    
    topo_out = cell(length(time_all),1);    % initialize output data
%     disp("# of time points: "+length(time_all));

    %% load boundary and generate inner boundary
    Coord = chan_coord_xy(:,1:2);
    load('borderCoords.mat','borderCoords');
    
    L = length(Coord); 
    inner_bound = generate_lining(borderCoords, 0.05);
    x2 = inner_bound(2,:);
    y2 = inner_bound(1,:);
    L2 = length(inner_bound);

    %normalize the data
    bord_maxX = max(borderCoords(2,:));
    bord_minX = min(borderCoords(2,:));
    bord_maxY = max(borderCoords(1,:));
    bord_minY = min(borderCoords(1,:));
    
    shift = [0 0 0 0];
    maxXnew = bord_maxX + shift(1);
    minXnew = bord_minX + shift(2);
    maxYnew = bord_maxY + shift(3);
    minYnew = bord_minY + shift(4);

    maxX = max(Coord(:,1));
    minX = min(Coord(:,1));
    maxY = max(Coord(:,2));
    minY = min(Coord(:,2));

    Coord(:,1) = ((maxXnew - minXnew) / (maxX - minX)) * (Coord(:,1) - maxX) + maxXnew;
    Coord(:,2) = ((maxYnew - minYnew) / (maxY - minY)) * (Coord(:,2) - maxY) + maxYnew;
    
    %% mesh the space
    X = Coord(:,1); 
    Y = Coord(:,2); 
%     border = 0; scatter = 0; edge = 1; extrapolateFlag = 0;
    % expand topoplot to border
    X(L+1:L+length(borderCoords)) = borderCoords(2,:);
    Y(L+1:L+length(borderCoords)) = borderCoords(1,:);
    
    %dummy
    x2(L2+1:L2+length(borderCoords)) = borderCoords(2,:);
    y2(L2+1:L2+length(borderCoords)) = borderCoords(1,:);
    
    resolution = 200;
    Xlin = linspace(min(X), max(X), resolution); 
    Ylin = linspace(min(Y), max(Y), resolution); 
    [xx,yy] = meshgrid(Xlin, Ylin);     
    
    %dummy plot
    Xlin2 = linspace(min(x2), max(x2), resolution); 
    Ylin2 = linspace(min(y2), max(y2), resolution); 
    [xx2,yy2] = meshgrid(Xlin2, Ylin2);
    
    z2 = ones(L2,1);

    % (don't extrapolate)
    z2(L2+1:L2+length(borderCoords)) = 100;
    
    %% loop through each time point
    for tt = 1:dt:round(size(rel_p,1)./tend)
        if mod(tt,300)==0; disp(tt); end
        
        Z = rel_p(tt,:)';
        
        % (don't extrapolate)
        Z(L+1:L+length(borderCoords)) = min(Z) + (max(Z) - min(Z)) / 2;

        z = griddata(X,Y,Z, xx,yy,'cubic'); 

        %compare dummy grid and actual grid to clean up edges
        % (edge==1) dummy grid for triming edges
        z0 = griddata(x2',y2', z2, xx2, yy2, 'cubic');
        for i = 1:resolution
            for j = 1:resolution
                if(z0(i,j) >= 95)
                    z(i,j) = NaN;
                end
            end
        end

        % clean up stray elements
        % (edge==1) dummy grid for triming edges
        for i = 2:resolution-1
            for j = 2:resolution-1
                if(~isnan(z(i,j)) && isnan(z(i+1,j)) && isnan(z(i-1,j)))
                    z(i,j) = NaN;
                end
            end
        end

        topo_out{tt} = reduce_resolution_circle_multithread(z,smooth,0.2); 

        % [~, MSGID] = lastwarn(); warning('off', MSGID);
    end

    TEND = toc(TSTART);
    disp(" > Elapsed time for this run: "+TEND/60+" mintues ("+TEND/3600+" hours).");
end