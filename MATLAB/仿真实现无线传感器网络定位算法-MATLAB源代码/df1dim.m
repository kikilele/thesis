function [df1] = df1dim(x)   
% Must accompany dlinmin.

% Defined in dlinmin.
global pcom
global xicom;
global nrdfun;

df  = feval(nrdfun, pcom + x.*xicom);
df1 = sum(df.*xicom);