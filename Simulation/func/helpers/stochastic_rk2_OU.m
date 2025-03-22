%
% Code for the solution of the
% differential equation:
%
%   dy/dt = f(t,y)
%


function [T, y] = stochastic_rk2_OU(f, init, Timepoints, numT, theta, mu, sigma)
    %  inputs
    y(:,1) = init;
    numT = numT - 1;
    h = (max(Timepoints) - min(Timepoints))/numT;
    
    T = min(Timepoints):h:max(Timepoints);
    dW = zeros(1,numT);
    %
    %  rk2 loop
    %
    for t = 1:numT
       Delta_Wk = sqrt(h)*randn;
       dW(1,t) = Delta_Wk;
       Sk = 1;
       if randn>0
          Sk = -1;
       end
       OU_k1 = h * theta * (mu - y(end-1,t)) + (Delta_Wk - Sk*sqrt(h))*sigma;
       SL_k1 = h * f( T(t), y(1:end,t));
       SL_k1(end-1,1) = OU_k1;
    
       OU_k2 = h * theta * (mu - y(end-1,t)) + (Delta_Wk + Sk*sqrt(h))*sigma;
       SL_k2 = h * f( T(t)+h, y(1:end,t) + SL_k1(1:end,1));
       SL_k2(end-1,1) = OU_k2;
    
       y(:,t+1) = y(:,t) + (SL_k1 + SL_k2)/2;
       y(end,t+1) = y(end,1)*exp(-theta*T(t)) + mu*(1-exp(-theta*T(t))) + sigma*sum(exp(-theta*(T(t)-T(1:t))).*dW(1:t));
    end
    
    y = y';
    T = T';
end