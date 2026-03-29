%  Used in the preparsing of for loops in the presence of multiple indexes
%  example
%  @for (indx,value)=enumerate({'S','I','R'})
%  	BIG_U@{indx} "welfare @{value}", C@{indx} "consumption (@{value})"
%  @end
% 
%  produces
% 
%  BIG_U1 "welfare S", C1 "consumption (S)"
%  BIG_U2 "welfare I", C2 "consumption (I)"
%  BIG_U3 "welfare R", C3 "consumption (R)"
%
%    Other uses of enumerate
%
%       rprt/enumerate
%