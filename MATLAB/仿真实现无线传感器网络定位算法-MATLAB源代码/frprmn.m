%|  FUNCTION:   frprmn
%|
%|  PURPOSE:  Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere 
%|            minimization is performed on a function func, using its gradient 
%|            as calculated by a routine dfunc. The convergence tolerance on 
%|            the function value is input as ftol. Returned quantities are p 
%|            (the location of the minimum), iter (the number of iterations 
%|            that were performed), and fret (the minimum value of the
%|            function). The routine linmin is called to perform line 
%|            minimizations.


function [p, iter, fret] = frprmn(p, ftol, func, dfunc, displayOption)

if (~exist('displayOption')),
    displayOption == 1;
end

%| 1. Constant definitions
ITMIN = 10;        % the minimum number of iterations before exit
ITMAX = 300;       % ITMAX is the maximum allowed number of iterations
EPS   = 1.0e-10;   % EPS is a small number to rectify the special case 
                   % of converging to exactly zero function value
if (displayOption == 1),
    disp('iter     mean xi  mini     p         xi   maxi     p         xi       fret');
    disp('----  ----------  ---------------------   ---------------------   --------');
end
                   
%| 2. Initializations
%disp('Initial:')
fp = feval(func, p);
xi = feval(dfunc, p);
exitCondition = 0;

g  = -xi;
h  = g;
xi = g;
extraLoops = 0;
maxExtraLoops = 10;

%| 3. Loop over iterations of minimization
for its=1:ITMAX,
   iter = its;
   [p, xi, fret] = dlinmin(p, xi, func, dfunc);
   absxi        = abs(xi);
   [minv, mini] = min(absxi);
   [maxv, maxi] = max(absxi);
   meanv        = mean(absxi);
   if displayOption==1,
       outStr = sprintf('%4d  %10.4g  %2d %7.4g %10.4g   %2d %7.4g %10.4g    %8.4g', ...
          iter, meanv, mini, p(mini), minv, maxi, p(maxi), maxv, fret);
       disp(outStr);
   end
     
   %| 4.  Normal exit condition
   if  ( 2.0*abs(fret-fp) <= ftol*( abs(fret) + abs(fp) + EPS )),
      if (its > ITMIN),
         exitCondition = 1;
         if (displayOption == 1),
             disp('Normal exit from frprmn.');
         end
         break;
      end
   end
   
   fp = fret;
   xi = feval(dfunc, p);
   gg = sum(g.^2);
   dgg = sum( (xi + g).*xi );   % This statement for Polak-Ribiere
   % dgg = sum( xi.^2);         % This statement for Fletcher-Reeves
   if gg == 0,            % Unlikely.  If gradient is exactly zero then 
      exitCondition = 2;   % we are already done.
      break;
   end
   gam = dgg/gg;
   g = -xi;
   h = g + gam.*h;
   xi = h;
end
if exitCondition == 0,
   disp('Too many iterations in frprmn');
end

