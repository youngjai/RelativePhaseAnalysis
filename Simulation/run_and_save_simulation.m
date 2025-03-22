% Run simulation with Kuramoto oscialltors and save output
% by Jae Hyung Woo; March 2025
addpath(genpath(pwd));
TSTART = tic;

model_type = "Kuramoto";

%%% Select stochastic model
%   OU: Ornstein-Uhlenbeck process
%   SPB3: supercritical/supercrifical bifurcation w/ 3rd order term
%   SPB5: subcritical/supercrifical bifurcation w/ 5th order term
%       for SPB, additional parameter 'r' is specified
stoch_model = "OU";
% stoch_model = "SPB5";  r = -0.18;
% stoch_model = "SPB3";  r = 0.1; 

%%% set phase delay term (tau)
tau_delay   = (0.12)*pi;

disp("Model: "+stoch_model);
disp("phase shift = "+tau_delay/pi+"pi");
%% Specify model & simulation parameters
time_ms = 20;     % averge window w/o sliding
% time_ms = 100;

numWorkers = parcluster('local').NumWorkers; 
disp("# workers = "+numWorkers);

numSim = numWorkers;  % # of init. cond. per parameter set
last_50_flag = 0;     % true if computing topo vector from last 50% of the data only

dt      = 0.001;  % time step for equilibration and simulation.
maxT1   = 30;      % length of time for equilibration with constant C0
maxT2   = 60;      % length of time for simulation with varying C0

S_set        = [1.8];   % global coupling strength
theta_set    = [6];   
sigma_set    = [1.5];   

mu_set       = sin(tau_delay);

disp("S set: "); disp(S_set);
disp("theta set: "); disp(theta_set);
disp("sigma set: "); disp(sigma_set);
disp("mu set: "); disp(mu_set);

num_of_configurations = length(theta_set)*length(sigma_set)*length(S_set)*length(mu_set);
disp("# of params set = "+num_of_configurations);
ParamsList = nan(num_of_configurations,4);
nCount = 0;
for q1 = 1:length(S_set)
    for q4 = 1:length(mu_set)
        for q2 = 1:length(theta_set)
            for q3 = 1:length(sigma_set)
                nCount = nCount + 1;
                ParamsList(nCount,:) = [S_set(q1), theta_set(q2), sigma_set(q3), mu_set(q4)];
            end
        end
    end
end

disp(ParamsList);
%% network coordinates
net = 'Lau114';    

switch net
    case 'Lau114'
        load('mat/Lau_connectivity_114.mat','MAT','coords');
        chan_coord_xy = [coords(:,2), -coords(:,1)];
end
N = length(MAT);

% load centroids (for regression on topoplot)
load('mat/combined_centroids_20240311.mat', 'C_comb', 'topo_idx_comb');
mask_topo = permute(C_comb.EOEC, [2,3,1]);
degree_num = sum(MAT,1);

%% Run simulations in parallel
parfevalOnAll(@warning,0,'off','all');

for ii = 1:size(ParamsList,1)
    iSTART = tic;
    disp("=====================================");
    disp("Set "+ii); disp(ParamsList(ii,:));
    if any(isnan(ParamsList(ii,:))); continue; end
    param = struct;
    
    % simulation param
    param.N       = N;      % number of oscillators
    param.dt      = dt;     % time step
    param.omega   = 10;     % intrinsic frequency, in Hz
    param.maxT1   = maxT1;  % equilibration time
    param.maxT2   = maxT2;  % simulation time
    param.noise   = 0;

    % model parm
    param.tau = tau_delay;    % phase shift
    param.S = ParamsList(ii,1);
    param.S_normN = 0;      % normalize S by N? Set zero as default
    param.W_std = 0;        % default 1/3; set to zero to fix all oscillators at omega
    
    % params for random process
    param.mu = ParamsList(ii,4);
    param.theta = ParamsList(ii,2);
    param.sigma = ParamsList(ii,3); 
    
    param.time_window = time_ms;    % ms
    param.time_moving = time_ms;
    param.last_50_flag = last_50_flag; % set to 1 if taking last 50% of the data only (for computing topo)
    disp("Time window: "+param.time_window+"ms");
    d_dir = "output/"+stoch_model+"_"+net+"/Kuramoto_"+time_ms+"ms/beta_"+num2str(tau_delay/pi,3)+"pi";
    if ~exist(d_dir,'dir'); mkdir(d_dir); end
    if contains(stoch_model,"SPB")
        disp("parameter r = "+r);
        param.r = r;
        dname = d_dir +"/"+net+"_DwellTime5_"+param.time_window+"ms_S"+param.S+"_norm"+param.S_normN+"_Wstd"+param.W_std+"_mu"+num2str(param.mu,4)+"_r"+r+"_Theta"+param.theta+"_Sigma"+param.sigma+"_n"+numSim;
    else
        dname = d_dir +"/"+net+"_DwellTime5_"+param.time_window+"ms_S"+param.S+"_norm"+param.S_normN+"_Wstd"+param.W_std+"_mu"+num2str(param.mu,4)+"_Theta"+param.theta+"_Sigma"+param.sigma+"_n"+numSim;
    end
    if param.last_50_flag; dname = dname + "_last50topo"; end
    dname = dname + ".mat";
    disp(dname);
    if exist(dname,'file')
        disp("File exists, skipping simulation: "+dname);
        continue;
    end
    
    %%% initial conditions w/ rseed
    clear OutDat
    parfor SEED = 1:numSim
        disp(" >> Simulating init. cond. "+SEED);       
        [X_prod, initC] = simulate_with_equil_Kuramoto(param, MAT, stoch_model, SEED);
        [rho1, pval1] = corr(degree_num', initC.W, 'type','Pearson');
        [rho2, pval2] = corr(degree_num', initC.W, 'type','Spearman');
        initC.Pearson_r_p = [rho1, pval1];  % correlation b/w degree # and intrinsic freq.
        initC.Spearman_r_p = [rho2, pval2];

        %%% calculate order parameter & Delta
        w = param.omega*2*pi;   % intrinsic freq.

        phase = X_prod(:, 1:N);
        deg2 = 1/N*sum(exp(1i*phase),2); % without amplitude term
        R_ord = abs(deg2);      % order parameter R
        delta = mod(phase(2:end,:) - phase(1:end-1,:)+pi,2*pi) - pi;
        Delta = w - mean(delta/dt,2);

        movwin = param.time_window*.001/param.dt;
        slidelen = param.time_moving*.001/param.dt;

        R_ord = moving_time_window(R_ord, slidelen, movwin); 
        Delta = moving_time_window(Delta, slidelen, movwin);
        D0    = moving_time_window(X_prod(:,end-1), slidelen, movwin);
        % subplot(3,1,1), plot(R2_ord), subplot(3,1,2), plot(Delta), subplot(3,1,3), plot(D0)
        
        % main processing func:
        disp("("+SEED+") Computing topo vector...");

        %%% compute topology matrices at each time point
        topo = compute_topo_func_Kuramoto(X_prod, chan_coord_xy, param, SEED);

        %%% K-mean clustering
        disp("("+SEED+") Calculating dwell time");

        %%% Regression indexing
        % edited by Youngjai Park, 2024.06.24. (Fri.) --------------------
        mask = permute(mask_topo,[3,1,2]);
        mask = mask(:,topo_idx_comb);
        K = size(mask,1);
        T_topo = size(topo,1);
        topo_vector = zeros([T_topo,length(topo_idx_comb)]);
        for t = 1:T_topo
            topo_vector(t,:) = topo{t}(topo_idx_comb);
        end

        [beta, epsilon] = cal_regression(mask, topo_vector);
        epsilon = mean(abs(epsilon), 1, 'omitnan'); % mean absolute error
        unit = param.time_window*.001;   % unit in sec
        [IDX, w_prop4] = cal_regression_clustering(beta, K, unit, 0);
        IDX4 = IDX; % 4 modes
        prop4 = cal_transition_prop_v2(IDX, K, unit);
        
        % infer Mode5 from epsilon values:
        error_prop = 0.15;
        ep_thresh = prctile(epsilon, (1 - error_prop)*100);
        IDX(epsilon>ep_thresh) = 5;
        [~, w_prop5] = cal_regression_clustering(beta, K+1, unit, 1000, IDX);
        prop5 = cal_transition_prop_v2(IDX, K+1, unit);

        OutDat(SEED) = prop5;
        OutDat(SEED).w_prop = w_prop5;  % weighted properties
        OutDat(SEED).IDX5 = IDX;    % 5 modes
        % OutDat(SEED).topo = topo;
        % ----------------------------------------------------------------
        % added beta & epsilon to OutDat, 2024.08.22
        OutDat(SEED).beta = beta;
        OutDat(SEED).epsilon = epsilon;
        OutDat(SEED).error_rate = error_prop;

        % Original info in 4 modes
        OutDat(SEED).Mode4_info = prop4;
        OutDat(SEED).Mode4_info.w_prop = w_prop4;   % weighted properties
        OutDat(SEED).Mode4_info.IDX = IDX4;   % 4 modes

        % model info
        OutDat(SEED).Model_info.R2 = R_ord;     % order parameter
        OutDat(SEED).Model_info.Delta = Delta;
        OutDat(SEED).Model_info.D0 = D0;
        OutDat(SEED).Model_info.initCond = initC; % initial condition

        disp(" >> "+SEED+": completed"); fprintf('\n');
    end    
    iEND = toc(iSTART);
    disp(" >> Elapsed time for this set: "+iEND/60+" minutes ("+iEND/3600+" hours).");
     
    save(dname, 'OutDat', 'param'); disp("File saved: "+dname);
    disp(datetime)
end
TEND = toc(TSTART);
disp(" >>>> Total elapsed time is "+TEND/60+" minutes ("+TEND/3600+" hours).");
