%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mod16qam.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = mod16qam(X,Fd,Fs,Eb);
% 16-QAM modulation
%----------------------
% Y = mod16qam(X,Fd,Fs,Eb)
%
% X - data stream
% Fd - data sampling rate
% Fs - Modulation signal samling rate (Fs/Fd integer)
% Eb - Average bit energy
%
M = 16;
half_d = sqrt(sqrt(0.4*Eb));
Y = half_d*dmodce(X,Fd,Fs,'qask',M)'; % QAM modulation

%qammod 
