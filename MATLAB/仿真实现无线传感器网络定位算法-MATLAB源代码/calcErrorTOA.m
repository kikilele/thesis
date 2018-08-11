%|
%| FUNTION:  calcErrorTOA
%|
%| PURPOSE:  Calculate the sum squared error between a model and the data
%|   when measurements are time-of-arrival (with Gaussian errors)
%|
%| ASSUMPTIONS:  
%| 1. that the blind device coordinates are put in 
%|    vectors which look like this:
%|                 [x1, x2, ... xn, y1, y2, ..., yn];
%|    where devices 1 to n are blind devices.
%|
%| 2. Assumes the dhat matrix is lower triangular.
%|
%| AUTHOR:  Neal Patwari 
%|   http://www.engin.umich.edu/~npatwari/
%| 

function [sumError] = calcErrorTOA(guessBlindLocs, parameters)

global refDevices;    % number of reference devices
global blindDevices;  % number of blind devices
global totalDevices;  % the total number of devices
global linearRefLocs; % actual coordinates of the reference devices
global dhat;          % estimated range between each pair of devices.
global funcEvals;     % counter for number of function evaluations.
funcEvals = funcEvals + 1;

%| 1.  Use these for the calculation of the relative location error term
x = [linearRefLocs(1:refDevices), guessBlindLocs(1:blindDevices)];
y = [linearRefLocs(refDevices+1:2*refDevices), ...
      guessBlindLocs(blindDevices+1:2*blindDevices)];

%| 2.  For each blind device, add in the error that is calculated between
%|     itself and all other blind or reference device.
sumError = 0;
for i = (refDevices+1) : totalDevices,
    j = 1:i-1;
    sumError = sumError + sum((sqrt((x(i)-x(j)).^2 + (y(i)-y(j)).^2) - dhat(i,j)).^2);
end

