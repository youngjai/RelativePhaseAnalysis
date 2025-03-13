function [coh, mpc, imc, pli, dpli, wpli, dwpli]=all_phase_measures_ht(bdata)
% this code utilize Hilbert transform
% input:
%   data: multivarite signal (time by channel)
%   segsize: segment size in bin
%   segmove: segment shift in bin
%   fs: sampling frequency
%
% output: matrix of Ch by Ch
%   coh: coherence
%   mpc: phase coherence
%   imc: imaginary coherence
%   pli: phase lag index
%   dpli: directed phase lag index


% Heonsoo Lee 2015.1.26
% updated 2015.8.30

[len, ch]=size(bdata);
a_sig=hilbert(bdata); % analytic signal
%% Measures using both amplitude and phase
sij=(a_sig'*a_sig)/len; % cross-spectra
sij=conj(sij);
coherency=sij./sqrt(diag(sij)*diag(sij)');

% 1. coherence
coh=abs(coherency);

% 3. imaginary coherence
imc=imag(coherency);

%% Measures only using phase

% phase=angle(a_sig); % extract phase
% a_sig2=exp(1i*phase); % reconstruct analytic signal after throwing out amplitude
a_sig2=a_sig./abs(a_sig);
sij2=(a_sig2'*a_sig2)/len; % cross-spectra
sij2=conj(sij2);

% 2. phase coherence
mpc=abs(sij2);

% 4. phase lag index & 5. weighted phase lag index
dpli=zeros(ch,ch);
for c1=1:ch-1
    for c2=c1+1:ch
        c_sig=a_sig(:,c1).*conj(a_sig(:,c2));
        ic_sig=imag(c_sig);
        dpli(c1,c2)=mean(sign(ic_sig));
        dpli(c2,c1)=-dpli(c1,c2);
        
        numer=mean(ic_sig);
        denom=mean(abs(ic_sig));
        dwpli(c1,c2)=numer./(denom+eps);
        dwpli(c2,c1)=-dwpli(c1,c2);
    end
end
pli=abs(dpli);
% pli(1:ch+1:end)=1;
wpli=abs(dwpli);






























