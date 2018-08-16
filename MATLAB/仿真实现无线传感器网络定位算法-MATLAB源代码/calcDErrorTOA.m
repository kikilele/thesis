%|
%| FUNTION:  calcDErrorTOA
%|
%| PURPOSE:  Calculate the derivative of the error between the model and the data
%|   when measurements are time-of-arrival (with Gaussian errors)
%|
%| ASSUMPTIONS:  
%| 1. that the blind (and reference) coordinates are put in a vector:
%|                 [x1, x2, ... xn, y1, y2, ..., yn];
%|
%| 2. Assumes the dhat matrix is lower triangular.
%| 
%| 3. the return value returns derivatives in this format:
%|       [d/dx1, d/dx2, ..., d/dxn,  d/dy1, d/dy2, ..., d/dyn];
%|
%| AUTHOR:  Neal Patwari 
%|   http://www.engin.umich.edu/~npatwari/
%| 

function [dError] = calcDErrorTOA(guessBlindLocs)

global refDevices;    % number of reference devices
global blindDevices;  % number of blind devices
global totalDevices;  % the total number of devices
global linearRefLocs; % actual coordinates of the reference devices
global dhat;          % estimated distance between devices based on the measured received power.
global dfuncEvals;    % counter for number of function evaluations.
dfuncEvals = dfuncEvals + 1;
TINY = 1e-8;

%| 1.  Use these for the calculation of the relative location error term
x = [linearRefLocs(1:refDevices), guessBlindLocs(1:blindDevices)];
y = [linearRefLocs(refDevices+1:2*refDevices), ...
     guessBlindLocs(blindDevices+1:2*blindDevices)];

%| 1.  Do the preliminary calculations here in order to save time in the
%|     next loop.
for k = (refDevices+1) : totalDevices,
   l = [1:k-1];
   modelDist       = max(TINY, sqrt(((x(k)-x(l)).^2 + (y(k)-y(l)).^2)));
   commonTerm(k,l) = 2.*(1 - dhat(k,l)./modelDist);
end
commonTerm(:,totalDevices) = zeros(totalDevices,1);

%| 2.  For each blind device, calculate the partial derivatives.
dFdx = zeros(1,totalDevices);
dFdy = zeros(1,totalDevices);
% row-wise error slopes
for k = (refDevices+1) : totalDevices,
    others   = 1:k-1;
    dFdx(k)  = sum(commonTerm(k,others).*(x(k)-x(others)));
    dFdy(k)  = sum(commonTerm(k,others).*(y(k)-y(others)));
end
% column-wise error slopes
for k = refDevices+1 : (totalDevices-1),
    others   = (k+1):totalDevices;
    dFdx(k)  = dFdx(k) + sum(commonTerm(others, k)'.*(x(k)-x(others)));
    dFdy(k)  = dFdy(k) + sum(commonTerm(others, k)'.*(y(k)-y(others)));
end

%| 3. Return the error in a vector.
dError = [dFdx((refDevices+1) : totalDevices), dFdy((refDevices+1) : totalDevices)];