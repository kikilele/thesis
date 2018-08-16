%|
%| FUNTION:  calcError
%|
%| PURPOSE:  Calculate the sum squared error between model and the data
%|   when measurements are RSS (with log-normal errors).
%|
%| ASSUMPTIONS:  
%| 1. that the blind (and the reference) device coordinates are put in 
%|    vectors which look like this:
%|                 [x1, x2, ... xn, y1, y2, ..., yn];
%|
%| 2. Assumes the dhat matrix is lower triangular.
%|
%| AUTHOR:  Neal Patwari 
%|   http://www.engin.umich.edu/~npatwari/
%| 

function [sumError] = calcError(guessBlindLocs)

global refDevices;    % number of reference devices
global blindDevices;  % number of blind devices
global totalDevices;  % the total number of devices
global linearRefLocs; % locations of the reference devices
global dhat;          % estimated distance between devices based on the measured
                      % received power.
global funcEvals;     % counter for number of function evaluations.
funcEvals = funcEvals + 1;
TINY      = 1e-5;

x = [linearRefLocs(1:refDevices), guessBlindLocs(1:blindDevices)];
y = [linearRefLocs(refDevices+1:2*refDevices), ...
      guessBlindLocs(blindDevices+1:2*blindDevices)];

%| 1.  For each blind device, add in the error that is calculated between
%|     itself and each other blind or reference device.
sumError = 0;
for i = refDevices+1 : totalDevices,
   j = 1:i-1;
   sumError = sumError + sum( log( max(TINY, ((x(i)-x(j)).^2 + (y(i)-y(j)).^2)) ...
      ./ (dhat(i,j).^2) ).^2);
end
