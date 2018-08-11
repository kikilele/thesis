%|
%| FUNTION:  calcDError
%|
%| PURPOSE:  Calculate the derivative of the error between the model and the data
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

function [dError] = calcDError(guessBlindLocs)

global refDevices;    % number of reference nodes
global blindDevices;  % number of blind devices
global totalDevices;  % the total number of devices
global linearRefLocs; % locations of the reference devices
global dhat;          % estimated distance between devices based on the measured
                      % received power.
global dfuncEvals;    % counter for number of function evaluations.
dfuncEvals = dfuncEvals + 1;
TINY       = 1e-5;

x = [linearRefLocs(1:refDevices), guessBlindLocs(1:blindDevices)];
y = [linearRefLocs(refDevices+1:2*refDevices), guessBlindLocs(blindDevices+1:2*blindDevices)];

%| 1.  Do the preliminary calculations here in order to save time in the
%|     next loop.
for k = refDevices+1 : totalDevices,
   l = [1:k-1];
   modelDistSqr    = max(TINY, (x(k)-x(l)).^2 + (y(k)-y(l)).^2);
   commonTerm(k,l) = log( modelDistSqr ./ (dhat(k,l).^2)) ./ modelDistSqr;
end

commonTerm(:,totalDevices) = zeros(totalDevices,1);

%| 2.  For each device, calculate the partial derivatives.
for k = refDevices+1 : totalDevices,
   dFdx(k)  = sum(commonTerm(k,1:k-1).*(x(k)-x(1:k-1))) + ...
      sum(commonTerm(k+1:totalDevices, k)'.*(x(k)-x(k+1:totalDevices)));
   dFdy(k)  = sum(commonTerm(k,1:k-1).*(y(k)-y(1:k-1))) + ...
      sum(commonTerm(k+1:totalDevices, k)'.*(y(k)-y(k+1:totalDevices)));
end

dError = [dFdx(refDevices+1:totalDevices), dFdy(refDevices+1:totalDevices)];
