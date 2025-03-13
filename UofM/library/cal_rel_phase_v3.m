% Calculate relative phase from time series
% Time length should be larger than channel length

% H = input time series (time x channel )
% rel_phase = relative phase original
% rel_phase_w = relative phase weighted
% rel_phase_t = relateve phase unweighted (triangle function)
% magni = complex magnitude
% R_theta = Global order parameter R
% T_theta = Global angle THETA
% m_theta = R exp(1i*THETA)

% 2023.03.29. Youngjai Park

function [rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta, e_theta] = cal_rel_phase_v3(H) 

        ch=size(H,2);  % channel size
        tp=size(H,1);  % time length

        z = hilbert(H);
%         To check the effect of the zero-padding
%         if tp < 200
%             z = hilbert(H,200);
%             z = z(1:tp,:);
%         else
%             z = hilbert(H);
%         end

        theta = angle(z); % The angles in theta are such that z = abs(z).*exp(i*theta).
        magni = abs(z);

        
        e_theta = exp(1i*theta); % exp(i*theta)
        m_theta = mean(e_theta,2);  %  m_theta = R exp(1i*THETA) = 1/N Sum{ exp(i*theta) } where N is number of channels
        T_theta = angle(m_theta);  % Global angle THETA
        R_theta = abs(m_theta);  % Global order parameter R

        T_theta_ch =  T_theta* ones(1,ch) ; % Prepare global angle THETA
        diff_global = theta- T_theta_ch ; % Calculate relative phase
        diff_global_rs =  mod( diff_global+pi,2*pi)-pi ; % Set range of relative phase into [-pi,pi] 

%         diff_global_rs2 = angle( exp(1i*diff_global) ); % Different expression to set range of relative phase into [-pi,pi]
%         % Compare two different expressions
%         figure(101);plot(diff_global_rs(:,1));hold on;plot(diff_global_rs2(:,1),'.'); % They are practically same 

        rel_phase = diff_global_rs;     % relative phase original version: discontinuity at pi and -pi
        rel_phase_w = sin(rel_phase);   % weighted relative phase: taking sine of the phase
        rel_phase_t = 2/pi* asin( sin ( rel_phase) ) ; % unweighted relative phase: triangle function transformation for relative phase 
%        [rel_phase_stan,rel_phase_stan_nm] = stan(rel_phase); 
        
end