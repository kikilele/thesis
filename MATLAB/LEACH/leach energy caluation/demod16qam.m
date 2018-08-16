%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demod16qam.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,Eb] = demod16qam(R,Fd,Fs);
% 16-QAM modulation
%----------------------
% [Y,Eb] = demod16qam(R,Fd,Fs)
%
% R - received signal (row vector)
% Fd - data sampling rate
% Fs - Modulation signal samling rate (Fs/Fd integer)
% Y - Output signal
% Eb - Average energy per received symbol
M = 16;
n = size(R,2);
Es = sum(abs(R).^2)/n;
Eb = Es/4;
half_d = sqrt(0.4*Eb);
% scatterplot(ynoisy,5,0,'b.'); % scatter plot of signal+noise
Y = ddemodce(R/half_d,Fd,Fs,'qask',M); % demodulated signal