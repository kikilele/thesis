%| FUNCTION: minBracket
%|
%| PURPOSE:  Given a function func, and given distinct initial 
%|   points ax and bx, this routine searches in the downhill 
%|   direction (defined by the function as evaluated at the 
%|   initial points) and returns new points ax, bx, cx that 
%|   bracket a minimum of the function. Also returned are the 
%|   function values at the three points, fa, fb, and fc.
%|
%| REFERENCE:  Numerical recipes in C "mnbrak"
%|
function [ax, bx, cx, fa, fb, fc] = minBracket(ax, bx, func)

%| 1. definitions
GOLD   = 1.618034;
GLIMIT = 100.0;
TINY   = 1.0e-20;

%|
fa   = feval(func, ax);
fb   = feval(func, bx);

%| 2. Switch roles for b and a so that we can go in a 
%|    downhill direction from b to a.
if fb > fa,
   temp = ax;   ax = bx;   bx = temp;
   temp = fb;   fb = fa;   fa = temp;
end
cx   = bx + GOLD*(bx-ax); % First guess for cx
fc   = feval(func, cx);

%| 3.  Keep returning here until we bracket
while fb > fc,
   r = (bx-ax) * (fb-fc);
   q = (bx-cx) * (fb-fa);
   u = bx - ((bx-cx)*q - (bx-ax)*r) / (2*nzSIGN(max(abs(q-r),TINY), q-r));
   
   %| 4. We won't go farther than ulim.  Test various possibilities.
   ulim = bx + GLIMIT*(cx-bx);
   if (bx-u)*(u-cx) > 0,  
      
      %| 4.1 parabolic u is between b and c: try it.
      fu = feval(func, u);
      if fu < fc,     % found a minimum between b and c
         ax = bx;  fa = fb;
         bx = u;   fb = fu;
         return;
      elseif fu > fb, % found a minimum between a and u
         cx = u;   fc = fu;
         return;
      end
           
      %| 4.2. parabolic fit was no use.  Use default magnification.
      u   = cx + GOLD*(cx-bx);
      fu  = feval(func, u);
      
   elseif (cx-u)*(u-ulim) > 0,
      %| 4.3.  Parabolic fit is between c and its allowed limit.
      fu = feval(func, u);
      if fu < fc,
         bx = cx;
         cx = u;
         u  = cx + GOLD*(cx-bx);
         fb = fc;
         fc = fu;
         fu = feval(func, fu);
      end
      
   elseif (u-ulim)*(ulim-cx) >= 0,
      %| 4.4. Limit parabolic u to its maximum allowed value.
      u = ulim;
      fu = feval(func, u);
   else
      %| 4.5 Reject parabolic u, use default magnification.
      u = cx + GOLD*(cx-bx);
      fu = feval(func, u);
   end
   
   %| 5. Eliminate oldest value and continue
   ax = bx; bx = cx; cx = u;
   fa = fb; fb = fc; fc = fu;
end
