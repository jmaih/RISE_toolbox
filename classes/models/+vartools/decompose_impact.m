function [omg,sig]=decompose_impact(C)
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


big=max(abs(C),[],1);
sig=diag(big);
omg=C./repmat(big,size(C,1),1);