%| nzSIGN  returns the abs() of the first argument if the 2nd argument is
%|    is >= 0.  Otherwise returnd the -abs() of the first argument.
%| 
%| PURPOSE:  written to match "numerical recipes in C" function SIGN.
%|    Unlike matlab's 'sign', a 0.0 value does not force the result to 0.0
%|

function [rval] = nzSIGN(valueArgument, signArgument)

if signArgument >= 0, 
   rval = abs(valueArgument);
else
   rval = -abs(valueArgument);
end
