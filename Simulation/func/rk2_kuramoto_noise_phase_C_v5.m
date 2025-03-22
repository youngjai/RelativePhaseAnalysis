% output:
% T is time, X is the resulting phase, XS is the resulting signal, W is the natural frequencies of the oscillators
% input:
% MAT is the connection matrix: default is 0
% noise is the noise: default is a gaussian distribution (normal distribution) with mean 0.1 and standard deviation 1
% phase is the extra phase term in the coupling: default is 0
% W is the natural frequency (in Hertz) of the oscillators: default is a gaussian distribution with mean 10 Hz and standard deviation 1  
% T is time period: default is 0:0.001:10 in second
% initV is the initial value: default is random

% example : 
% load(connection_matrix_gong_MAT.mat);
% [T,X] = rk2_kuramoto_noise_phase(MAT);

% 2014.3.6. written by Joon-Young Moon
% 2019.2.22 modified by Joon-Young Moon
% 2024.10.11 modified by Jae Hyung Woo
% This code requires rk2.m

function [T, X] = rk2_kuramoto_noise_phase_C_v5(MAT, Beta, C0, W, noise, T, initV)
%     tic
    N = length(MAT); % number of ocillators 
    W_rad = W.*2*pi;
    numT = length(T);
    
    [m_x, m_y, m_z] = find(MAT);
    sumMAT = sum(MAT);  % sum of columns (node degrees)
      
    [T, X] = rk2(@(T,X) J(X, W_rad, Beta, C0, N, sumMAT, m_x, m_y, m_z, noise), initV, T, numT);   
%     ET = toc;   disp(ET/60+" minutes");
    return
end

function dx = J(X, W, Beta, C0, N, sumMAT, m_x, m_y, m_z, noise)
         dx = zeros(N,1);
         XT_list = m_z .* sin( X(m_x) - X(m_y) - Beta ); 
         term1 = sum(sparse(m_x, m_y, XT_list, N, N));
         term2 = sumMAT.*C0;

         dx(:) = W(:) + term1' + noise.*randn(N,1) + term2' ;
    return
end
