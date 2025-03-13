function [coh, pc, imc, pli, dpli, psi, wpli, dwpli, f]=all_phase_measures_ft(data, segsize, segmove, fs, df)
% this code utilize Fourier transform 
% input:
%   data: multivarite signal (time by channel)
%   segsize: segment size in bin
%   segmove: segment shift in bin
%   fs: sampling frequency
%
% output:  matrix of Ch by Ch by frequency
%   coh: coherence
%   pc: phase coherence
%   imc: imaginary coherence
%   pli: phase lag index
%   dpli: directed phase lag index
%   psi: phase slope index


% Heonsoo Lee 2015.1.26

ch=size(data,2); % number of channels
f=linspace(0, fs/2, segsize/2+1);

segnum=floor((length(data)-segsize)/segmove+1); % number of segments
hw=hanning(segsize);
hw=repmat(hw,1,ch); % hamming window

for u=1:segnum % for each segment
    segdata=data(segmove*(u-1)+1:segmove*(u-1)+segsize,:); % segment data
    segdata_hw=detrend(segdata).*hw; % hamming windowed segment data
    y=fft(segdata_hw); % fourier spectra
    y=y(1:end/2+1,:); % take half
    for ff=1:size(y,1) % for each frequency
        yij(:,:,ff,u)=conj(y(ff,:)'*y(ff,:)); % cross-spectra
    end
end

%%% coherency %%%
sij=mean(yij,4); % average over segments (time)
for ff=1:size(sij,3)
    coherency(:,:,ff)=sij(:,:,ff)./sqrt(diag(sij(:,:,ff))*diag(sij(:,:,ff))');
end

% 1. coh
coh=abs(coherency);

% 2. pc
% phi=angle(yij);yij2=exp(1i*phi);
yij2=yij./abs(yij);
pc=abs(mean(yij2,4));

% 3. imc
imc=imag(coherency);

% 4. phase lag index & 5. directed phase lag index
dpli=mean(sign(imag(yij)),4);
pli=abs(dpli);
% pli(1:ch+1:end)=1;

% 6. phase slope index
if nargin <5
    df=1;
end
psi=imag(conj(coherency(:,:,1:end-df)).*coherency(:,:,1+df:end));

% 7. directed weighted phase lag index
numer=mean(imag(yij),4); % average of imaginary
denom=mean(abs(imag(yij)),4); % average of abs of imaginary
dwpli=numer./(denom+eps);
wpli=abs(dwpli);

































