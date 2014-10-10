function [r,Jac,retcode]=ss_residuals(ss_i,resid_func,func_jac,x_ss,pp_i,def_i)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


retcode=0;
Jac=[];

y_=ss_i;
% input list is always 'y'    'x'    'ss'    'param'    'sparam'    'def'    's0'    's1'
r=utils.code.evaluate_functions(resid_func,y_,x_ss,ss_i,pp_i,[],def_i,[],[]);
if nargout>1
    Jac=func_jac(y_,x_ss,ss_i,pp_i,[],def_i,[],[]);
end

if retcode && obj(1).options.debug
    utils.error.decipher(retcode)
end
end