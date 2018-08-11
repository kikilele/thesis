%| FUNCTION: dlinmin
%|
%| PURPOSE:  Given an n-dimensional point p[1..n] and an 
%|   n-dimensional direction xi[1..n], moves and resets p to where 
%|   the function func(p) takes on a minimum along the direction xi 
%|   from p, and replaces xi by the actual vector displacement that 
%|   p was moved. Also returns as fret the value of func at the returned 
%|   location p. This is actually all accomplished by calling the
%|   routines mnbrak and brent.
%|
%| REFERENCE:  Numerical recipes in C
%|

function [p, xi, fret] = dlinmin(p, xi, func, dfunc)

TOL    = 2.0e-4;  % Tolerance passed to brent.

global pcom xicom nrfunc nrdfun;
nrfunc = func;
nrdfun = dfunc;
pcom   = p;
xicom  = xi;

ax     = 0.0;  % Initial guess for brackets.
xx     = 2.0;
[ax, xx, bx] = minBracket(ax, xx, 'f1dim');
[fret, xmin] = dbrent(ax,xx,bx,'f1dim', 'df1dim',TOL);

xi     = xi.*xmin;
p      = p + xi;
