clear all
clc
%% create function
func=@(a,b,c,d)[exp(a+2*log(b+c)-a*atan(b*c));exp(-a*atan(b*c));d*exp(a+2*log(b+c))];
%% initialize arguments
a=rise_sad('a');b=rise_sad('b');c=rise_sad('c');d=rise_sad('d');
%% create tree
tree_=func(a,b,c,d);
%% compute derivatives... and register the number of calls to each node in the tree 
[dd,references]=diff(tree_,{a,b,c,d});
%%
[Jac,dd,references]=rise_sad.jacobian(func,{a,b,c,d},{a,b})
%% Jacobian with strings
[Jac,dd,references]=rise_sad.jacobian({'exp(a+2*log(b+c)-a*atan(b*c))'
    'exp(-a*atan(b*c))'
    'd*exp(a+2*log(b+c))'},{'a','b','c','d'},{'a','b'})
%% Jacobian with strings and no wrt
profile off
profile on
[Jac,dd,references]=rise_sad.jacobian({'exp(a+2*log(b+c)-a*atan(b*c))'
    'exp(-a*atan(b*c))'
    'd*exp(a+2*log(b+c))'},{'a','b','c','d'})
profile off
profile viewer
%% Jacobian with tree
[Jac,dd,references]=rise_sad.jacobian(tree_,{'a','b','c','d'})
%% Jacobian with tree and no arguments must have wrt in that case
[Jac,dd,references]=rise_sad.jacobian(tree_,[],{'a','b'})
%% Hessian
[H,Jac_char,references]=rise_sad.hessian('exp(a+2*log(b+c)-a*atan(b*c))',{'a','b','c','d'},{'a','b','c','d'})