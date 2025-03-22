% function for running simulation with specified parameters
% output is matrix of node phases
function [X_prod, initC] = simulate_with_equil_Kuramoto(param, MAT, stoch_model, SEED)
    if ~isfield(param,'S_normN')
        param.S_normN = 0;
    end
    N       = param.N;      % number of oscillators (nodes)
    dt      = param.dt;     % time step
    maxT1   = param.maxT1;  % simulation time for equilibration process
    maxT2   = param.maxT2;  % simulation time for stochastic process for d_0
    
    Tau = param.tau;
    noise = param.noise;
    Theta = param.theta;    % drift rate
    Sigma = param.sigma;    % std. of random walk
    Mu = param.mu;          % baseline (+ initial) d_0 value
    S = param.S;            % global coupling strength
    
    if param.S_normN
        S_MAT = S.*MAT / N;  % if normalizng coupling strength by N
    else
        S_MAT = S.*MAT;      % no normalization; this is default
    end

    rng(SEED, 'threefry');
    % set initial values around a unit circle:
    init_theta = mod(pi + 0.1*pi*randn(1,N), 2*pi);     % Gaussin with std=0.1pi
    
    initC = struct;     % initial condition
    initC.init_theta = init_theta';

    C0_init = Mu;
    W = (param.omega).*ones(1,N) + param.W_std * randn(1,N);   % in Hz (will be multiplied by 2pi in below func)
    initC.W = W';                                     % % 99.7% (3std) falls between [9, 11]
    T1 = 0:dt:maxT1;
    T2 = 0:dt:maxT2;

    % Equilibration: no random process 
    initV = init_theta;
    [~, X_equil] = rk2_kuramoto_noise_phase_C_v5(S_MAT, Tau, C0_init, W, noise, T1, initV);  
    
    % Reshape initV to include D0.
    initV = zeros(1, N+2);
    initV(1:N) = X_equil(end,:);
    initV(end-1) = C0_init;     % Last indices of initV are simulation result of D0.
    initV(end) = C0_init;       % Last indices of initV are analytical result of D0; not used

    % Simulation with v0
    switch stoch_model
        case 'OU'
            [~, X_prod] = rk2_OU_Kuramoto_network_phase_list_C0(S_MAT, Tau, noise, W, T2, initV, Theta, Mu, Sigma);
        case 'SPB3'
            Order = 3;
            [~, X_prod] = rk2_SPB_Kuramoto_network_phase_list_C0(S_MAT, Tau, noise, W, T2, initV, Theta, Mu, Sigma, param.r, Order);
        case 'SPB5'
            Order = 5;
            [~, X_prod] = rk2_SPB_Kuramoto_network_phase_list_C0(S_MAT, Tau, noise, W, T2, initV, Theta, Mu, Sigma, param.r, Order);
    end
end