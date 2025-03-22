% 2024.10.11 Produced by J.H. Woo, freely distributable
%   simulates stochastic process on C0 values with subcritical pitchfork bifurcation (SPB)
%     
% OUTPUT
%   T is time, X is the resulting signal
% INPUT
%   MAT is the connection matrix
%   T is time
%   initV is the initial condition
%   Tau is phase shift (delay) term
%   initV is the initial condition including phase and initial C0
%       C0 is external input term
%   theta is drift rate in SPB stochastic process
%   mu is fixed point in SPB stochastic process
%   sigma is std. of random noise in stochastic process
%   spb_r is the parameter r in SPB
%   Order: specify as either 3rd-order (3) or 5th-order (5)

%   Note: ratio of alpha frequency vs. delay phase is defined as:  Tau = 2pi*delay/100ms   
%   Beta = 0.2*pi; in case of delay of 10ms

function [T, X] = rk2_SPB_Kuramoto_network_phase_list_C0(MAT, Tau, noise, W, T, initV, theta, mu, sigma, spb_r, Order)
% tic
    N = length(MAT);    % number of oscillators
    W_rad = W.*2*pi;    % convert HZ into radians
    numT = length(T);   % number of time points

    [m_x, m_y, m_z] = find(MAT);    % indices for nonzero connection
    sumMAT = sum(MAT);  % node degrees

    switch Order
        case 3
            [T, X] = stochastic_rk2_SPB_3rd(@(T,X) J(X, W_rad, Tau, N, sumMAT, m_x, m_y, m_z, noise), initV, T, numT, theta, mu, sigma, spb_r);
        case 5
            [T, X] = stochastic_rk2_SPB_5th(@(T,X) J(X, W_rad, Tau, N, sumMAT, m_x, m_y, m_z, noise), initV, T, numT, theta, mu, sigma, spb_r);
    end
% toc
    return
end

function dx = J(X0, W, Tau, N, sumMAT, m_x, m_y, m_z, noise)
        dx = zeros(N+2, 1);
        
        X  = X0(1:N);
        C0 = X0(N+1);
        % C0_exact = X0(end);   % not used

        XT_list = m_z .* sin( X(m_x) - X(m_y) - Tau ); 
        term1 = sum( sparse(m_x, m_y, XT_list, N, N) );
        term2 = sumMAT.*C0;

        dx(1:N) = W(:) + term1' + noise.*randn(N,1) + term2' ;
    return
end
