%
% Code for the solution of the
% differential equation:
%
%   dy/dt = f(t,y)
%
function [T, y] = stochastic_rk2_SPB_3rd(f, init, Timepoints, numT, theta, mu, sigma, r)
    %  inputs
    y(:,1) = init;
    numT = numT - 1;
    h = ( max(Timepoints) - min(Timepoints) ) / numT;
    
    T = min(Timepoints) : h : max(Timepoints);
    dW = zeros(1,numT);
    
    %  rk2 loop
    for t = 1:numT
       Delta_Wk = sqrt(h)*randn;
       dW(1,t) = Delta_Wk;
       Sk = 1;
       if randn>0
          Sk = -1;  % randomly determine direction (sign) of the stochastic step
       end

       % 3-rd order SPB
       d0_k1 = h * theta*( r * (y(end-1,t) - mu) - (y(end-1,t) - mu)^3 ) + sigma * (Delta_Wk - Sk*sqrt(h));
       SL_k1 = h * f( T(t), y(1:end,t));
       SL_k1(end-1,1) = d0_k1;
       
       % 3-rd order SPB
       d0_k2 = h * theta*( r * (y(end-1,t) - mu) - (y(end-1,t) - mu)^3) + sigma * (Delta_Wk + Sk*sqrt(h));
       SL_k2 = h * f( T(t)+h, y(1:end,t) + SL_k1(1:end,1));
       SL_k2(end-1,1) = d0_k2;
       
       y(:,t+1) = y(:,t) + (SL_k1 + SL_k2)/2;
    end
    
    y = y';
    T = T';
end